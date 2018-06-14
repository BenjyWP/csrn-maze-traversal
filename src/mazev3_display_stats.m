function mazev3_display_stats(net, dW, dww, dWs, MAZES, EE,p)
%%% Covariance Matricies and their eigenvalues 
%%% should be able to shed light on what's going on

%% K is error conariance matrix
% e = eig(net.K);
% ee=sort(e, 'descend');
% plot(ee(16:end));
% pause;

%% C*K*C' + R is the measurement covariance matrix
% e= eig(net.C*net.K*net.C' + net.R);
% ee=sort(e, 'descend');
% plot(ee);
% pause;

%%Look at the "weighting factors"
B = (net.R * inv(net.C*net.K*net.C' + net.R));
% surf(B(:,:));
% pause(2);
return

%% Let's compare how the calculated changes in w affect the Error
TotdeltaW = net.C'*net.alpha; % using simple Backprop
deltaW  = -0.00000000001*reshape(TotdeltaW(1:net.n*(net.n+net.m-1)), net.n, (net.n+net.m-1));
deltaww = -0.00000000001*TotdeltaW(net.n*(net.n+net.m-1)+1:net.n*(net.n+net.m-1)+net.n);
deltaWs = -0.0000000001*TotdeltaW(end);

% make a new net and change weights
net1 = net;
net1.c_W(net1.m:net1.m+net1.n-1,:) = net1.c_W(net1.m:net1.m+net1.n-1,:) + deltaW;
net1.c_ww = net1.c_ww + deltaww;
net1.Ws = net1.Ws + deltaWs;
% calculate forward to get the new error
TOT_TEST_E= 0;
for test_maze_num=1:size(MAZES,2)
    maze=MAZES{test_maze_num};        
    net1 = mazev3_assign_external_inputs(maze, net1);
    %%% forward calculation
    [y, v, net1.store_y]=MAZE_CELLULAR_SRN(net1.c_W, net1.c_ww, net1.c_x, net1.n, net1.m-1, maze.size_x, p);
    [TEST_E , c_J, c_DeltaJ]= mazev3_read_output(maze, net1, p);
    TOT_TEST_E = TOT_TEST_E + TEST_E;
end    

