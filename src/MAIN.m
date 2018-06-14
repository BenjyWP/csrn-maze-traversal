%==========================================================================
% EKF solution w/ Added comments
%==========================================================================
%--------------------------------------------------------------------------
% Function:  Mazev3_main_M5ekf_000.m
% Description: This version of the mazev3 code has is basically a cleaned
%              up version of Roman's code.  The ALR option has been removed
%              and only the Kalman method(EKF) remains.  Additional
%              comments have been added.
% Purpose:  The purpose of this code version is to clean up to code,
%           removing portions that we are not using...and to add addition 
%           comments to facilitate a better understanding of the code.
% Noteable Features:
%       (1) uses noise annealing

% Output:  
%       (1) Plot: output & target mazes for test mazes
%       (2) Plot: SSE vs epochs
%       (3) Plot: Goodness of Solution vs epochs
%       (4) Plot: Measurement Noise vs epochs
%       (5) Plot: Selected Weights vs epochs
%--------------------------------------------------------------------------
% Basis: mazev3
% Operand: maze
% Operand size: 5x5
% Trianing Method: EKF
% Training: variable # of mazes
% Testing: variable # of mazes
% VerID:  001
% Date: 4/18/08
% -------------------------------------------------------------------------
% Version History:
%
% Basis:  This code is based on Roman Illiad's mazev3 code:  Roman's code
% uses the CSRN network to solve the basic maze problem.
%
% 000 - initial version
% -------------------------------------------------------------------------
% Modification History:
% 02/01/08:  original cleanup & comments
% 03/26/08:  added plots for target maze, SSEvsEpochs, GOSvsEphochs
%            and MNOISEvsEpochs.
% 04/04/08:  added plots for selected weights vs epochs
% 04/18/08:  added header info.
% -------------------------------------------------------------------------
% Primary Researcher: K.Iftekharrudin
% Lab: ISIP
% Authors:
% (1) Roman Illiad, original mazev3
% (2) Keith Anderson
%--------------------------------------------------------------------------

clear all;
% close all;
clc;

%kalman variables
%net.K = covariance matrix
%net.Q = process noise = 0 for this application
%net.R = measurement noise, annealed for this application
%net.c_W = weights
%net.c_ww = bias node weights
%net.c_Ws = scaling factor weight
%net.c_x = external inputs
%net.C = jacobian matrix

MAZE_SIZE = 5;      % maze size
NUM_MAZES = 5;      % number of mazes for training/testing  (Roman used 6)

%generate training & testing mazes
MAZES = mazev3_generate_random_mazes(NUM_MAZES, MAZE_SIZE, MAZE_SIZE, 'TRAINING MAZES');
TEST_MAZES = mazev3_generate_random_mazes(NUM_MAZES, MAZE_SIZE, MAZE_SIZE, 'TESTING MAZES');
%TEST_MAZES=MAZES;
%===== generate initial values for NN
%randomly generate mulitple sets of initial values for NN.  Test to see
%which set gives best result.  Return this set for our initial NN values.
NUM_RAND_NETS = 40;                     %how many rand nets to generate
[net, MIN_EE] = mazev3_generate_random_net(MAZES , NUM_RAND_NETS);

net.METHOD = 'KALMAN';                  % Learning method
BATCH_KALMAN = 0;                       % switch to batch Kalman

%initialize measurement noise & annealing function
maze=MAZES{1};
net.R = net.R*eye(size(MAZES,2)*maze.size_y*maze.size_x,size(MAZES,2)*maze.size_y*maze.size_x);
net.inline_r = inline('0.003*log(0.001*E+1)', 'E');  %annealing function

init_net = net;                         % save the initial net

%===== set up other params
p = net.MAX_STEPS;                      % initial "fast" interations
t_last_change=0;                        % counter used for incrementing p
override=1;                             % helps with stopping cond'n
EPOCHS = 10;                            % max epochs

ERRORS= [];                             % matrix to store errors and other info
MIN_TOT_TEST_E = Inf;
MAX_TOT_TEST_GOODNESS = 0;
best_t=0;
best_good_t = 0;
myWeights=zeros(7,EPOCHS);              %weight matrix to hold capture the weights as they change over epochs

tic
for t=1:EPOCHS
    %init some variables
    EE=0;                                       % overall sum-squared error(across all mazes)
    DJ = zeros(maze.size_y,maze.size_x);        % overall cell error(across all mazes)           
    TOT_PCT = 0;                                % percentage of correct directions(for all mazes)

    % Testing section -----------------------------------------------------
    % (Training????)
    for maze_num=1:size(MAZES,2)
        maze=MAZES{maze_num};                   %capture current maze
        net.t=t;                                %capture current epoch        
        
        %assign external node inputs(net.c_x)
        net = mazev3_assign_external_inputs(maze, net);     
        
        %----- forward calculation
        % y = output for each node(1 cell = 16 nodes). 16x49 cells = 784
        %     (last column of net.store_y)
        % v = output for each OUTPUT node. 5 x 49 cells = 245
        %      (have not verified yet)
        % net.store_y = output for each node over p iterations. 784 x 10
        [y, v, net.store_y]=MAZE_CELLULAR_SRN(net.c_W, net.c_ww, net.c_x, net.n, net.m-1, maze.size_x, p,0,0,0,0);
        
        %----- extract cell outputs from net.store_y & compute cell error
        % c_J = the output for the cells of the maze 
        % c_DeltaJ = the error for the current maze(ie the diff between
        %            maze output and the target maze.
        % E = sum-squared error for current maze
        [E , c_J, c_DeltaJ]=mazev3_read_output(maze, net, p);
             
        %----- compute statistics
        EE = EE+E;                                          %compute overall SSE(across all mazes)
        DJ = DJ + c_DeltaJ;                                 %compute cell error(across all mazes)    
        pct = mazev3_goodness_of_solution(maze.A, c_J);     %percentage of correct directions for this maze 
        TOT_PCT = TOT_PCT + pct;                            %percentage of correct directions(for all mazes)
        
        %----- reshape error array, c_DeltaJ
        %c_DeltaJ = F_J = maze(cell) error...ie c_J - Maze.A
        %reshape c_DeltaJ into a vector...[49x1]
        F_J = reshape( c_DeltaJ',maze.size_x*maze.size_y,1);    %non-zero only the first iteration of SRN
        
        %The backpropagation computation computes
        %        d(error)/dW = dE/dW. 
        %For our EKF algorithm we need the Jacobian, which is 
        %        d(output)/dw = dC_J/dW.
        %By forcing the error, c_DeltaJ = F_J, to 1, dE/dW becomes
        %          dE/dW = dC_J/dw.
        %This "trick" allows us to use the backprop code to calculate the
        %Jacabian, dC_J/dw, instead of dE/dW.
        if strcmp(net.METHOD , 'KALMAN')==1
            F_J(:,:)=1;     %this is for Kalman filter, so that C is the Jacobian, not F_W        
        end 
        
        [F_W, F_ww, F_Ws, F_x, F_Y, net.store_f_y] = MAZE_CELLULAR_SRN_FEEDBACK(net.c_W, net.c_ww,  net.Ws, net.n, net.m-1, maze.size_x, p, net.store_y, F_J);
        
        %--------------------------- KALMAN Method
        if strcmp(net.METHOD , 'KALMAN')==1            
            if maze_num == size(MAZES,2)
                ERRORS = [ERRORS; t EE net.R(1,1) 0 0 0];                                     %build up ERROR array
               
                net = mazev3_KalmanAddRowToJacobian(net, maze, c_DeltaJ, F_W, F_ww, F_Ws);    %form jacobian
             
                [TotalDW,net]=mazev3_KalmanStep(net, maze, DJ,EE,ERRORS, []);                 %perform EKF calc.                      
              
                dW  = reshape(TotalDW(1:net.n*(net.n+net.m-1)), net.n, (net.n+net.m-1));      %capture chang in weights
                dww = TotalDW(net.n*(net.n+net.m-1)+1:net.n*(net.n+net.m-1)+net.n);
                dWs = TotalDW(end);
                %mazev3_display_stats(net, dW, dww, dWs, MAZES, EE, p);
                
                %%% Different Weight Adjustment for C++ version
                net.c_W(net.m:net.m+net.n-1,:) = net.c_W(net.m:net.m+net.n-1,:) + dW;         %update weights
                net.c_ww = net.c_ww + dww;
                net.Ws = net.Ws + dWs;
                net.C = [];                             %clear jacobian matrix
                net.alpha=[];                           %clear inovation matrix
            else
                if (BATCH_KALMAN~=1)
                    %add row to net.C
                    net = mazev3_KalmanAddRowToJacobian(net, maze, c_DeltaJ, F_W, F_ww, F_Ws);    %form jacobian
                end
            end
        end
    end  %maze loop for training
   
    %----- capture weights for this epoch
    myWeights(1,t)= t;                                      %capture record epoch #
    myWeights(2:6,t)= net.c_ww';                            %capture recurrent weights 
    myWeights(7,t)= net.Ws;                                 %   and scaling weight
    
    
    %% Testing phase ------------------------------------------------------
    TOT_TEST_E= 0;
    TOT_TEST_PCT=0;
    
    for test_maze_num=1:size(TEST_MAZES,2)                          %----- maze loop
        maze=TEST_MAZES{test_maze_num};                             %capture current maze    
        net = mazev3_assign_external_inputs(maze, net);             %assign inputs
        %%% forward calculation
        [y, v, net.store_y]=MAZE_CELLULAR_SRN(net.c_W, net.c_ww, net.c_x, net.n, net.m-1, maze.size_x, p,0,0,0,0);
        % y: 1 by 735, v: 245 by 1, net.store_y: 735 by 18
        
        [TEST_E , c_J, c_DeltaJ]= mazev3_read_output(maze, net, p); %extract error data
        
        pct = mazev3_goodness_of_solution(maze.A, c_J);             %calc. stats
        TOT_TEST_PCT = TOT_TEST_PCT + pct;
        TOT_TEST_E = TOT_TEST_E + TEST_E;
    end %maze loop for testing
    
    %complete error info for this epoch
    ERRORS(end, 4) = TOT_TEST_E;                                    
    ERRORS(end, 5) = TOT_PCT/size(MAZES,2);
    ERRORS(end, 6) = TOT_TEST_PCT/size(TEST_MAZES,2);
    
    %Save the best net based on testing error 
    if MIN_TOT_TEST_E > TOT_TEST_E
        MIN_TOT_TEST_E = TOT_TEST_E;
        best_gene_net = net;        
        best_t=t;
    end

    %save the best net based on goodness of navigation in test mazes
    if MAX_TOT_TEST_GOODNESS < TOT_TEST_PCT
       MAX_TOT_TEST_GOODNESS =  TOT_TEST_PCT;
       best_good_net = net;           
       best_good_t=t;
    end
    
    %----- early stopping code -----
    %[t EE net.R(1,1) TOT_PCT/size(MAZES,2) TOT_TEST_PCT/size(TEST_MAZES,2) net.Lr_W]
    if override~=1        
        [r,override]=mazev3_stoppingcondition(net,maze, ERRORS);
        if  and(r==1,override~=1)
            ['Stopping...']
            break;
        end
    end
    if (TOT_PCT/size(MAZES,2))>80
        ['Stopping...']
        break;
    end
    fprintf('Epoch: %d\n', t);
end %t, epoch loop
toc

% Result Section ----------------------------------------------------------

%----- select best generalizing net
% for the best generalizing network, for each test maze, do the forward
% calculation, extract the output info, compute the error, and display the
% maze.
TOT_E=0;
for maze_num=1:size(TEST_MAZES,2)                                           %maze loop
    maze=TEST_MAZES{maze_num};                                              %capture current maze       
    best_gene_net = mazev3_assign_external_inputs(maze, best_gene_net);     %assign external inputs  
    %%% forward calculation
    [y, v, best_gene_net.store_y]=MAZE_CELLULAR_SRN(best_gene_net.c_W, best_gene_net.c_ww, best_gene_net.c_x, best_gene_net.n,...
        best_gene_net.m-1, maze.size_x, p,0,0,0,0);
    [E , c_J, c_DeltaJ]= mazev3_read_output(maze, best_gene_net, p);        %extract output data
    TOT_E=TOT_E+E;                                                          %compute error
    
    maze_v3_display_maze(c_J,1);                                            %display maze
    maze_v3_display_maze(maze.A,1);
end

%----- plot error vs epochs ----
figure;
hold on;
xlabel('epochs')
ylabel('SSE')
title('SSE - over all test mazes')
plot(ERRORS(:,1),ERRORS(:,4),'k');              %SSE over all mazes(TEST mazes)
hold off;

%----- plot goodness of solution vs epochs ----
% figure;
% hold on;
% xlabel('epochs')
% ylabel('%')
% title('Average Goodness of Solution')
% plot(ERRORS(:,1),ERRORS(:,6),'y');              %Ave Goodness of Solution(TEST mazes)
% hold off;

%----- plot Measurment Noise vs epochs ----
% figure;
% hold on;
% xlabel('epochs')
% ylabel('Co-variance')
% title('Measurement Noise Covariance')
% plot(ERRORS(:,1),ERRORS(:,3));              %Measurement Noise Covariance(TEST mazes)
% hold off;

figure;
hold on;
xlabel('epochs')
ylabel('scaling weight')
title('Weights across epochs')
%plot(myWeights(1,:),myWeights(7,:));                   %scaling weight
subplot(3,2,1); plot(myWeights(1,:),myWeights(2,:));    %bias weights
subplot(3,2,2); plot(myWeights(1,:),myWeights(3,:));
subplot(3,2,3); plot(myWeights(1,:),myWeights(4,:));
subplot(3,2,4); plot(myWeights(1,:),myWeights(5,:));
subplot(3,2,5); plot(myWeights(1,:),myWeights(6,:));
subplot(3,2,6); plot(myWeights(1,:),myWeights(7,:));    %scaling weight
hold off;