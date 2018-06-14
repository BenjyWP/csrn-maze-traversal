function [net] = mazev3_initialize_net(objMaze)
SPARSE = 1;
net.n=5;   % number of recurrent nodes - outputs 4
net.m=net.n + 4 + 3;  % number of net inputs + bias 12;  original 4+3
net.N=net.n + 4 + 3;  % number of hidden units + inputs _ bias 12;  original 4+3
net.Cell_Nodes = net.N+net.n; % total number of cell nodes
net.MAX_STEPS = 10;     %maximum recurrent steps
net.t = 0;  %current training epoch
%Learning methods
%net.METHOD = 'ALR'; 
net.METHOD = 'KALMAN';
% Generate Weight Matrices 
net.Ws=40;
%% Kalman Filter Parameters - exp 2
net.R = 10*eye(1,1);
net.K = 0.00001 .* eye(1 + net.n + net.n*(net.n+net.m-1),1 +  net.n + net.n*(net.n+net.m-1));
net.Q = 0.*0.00001 .* eye(1 + net.n + net.n*(net.n+net.m-1),1 + net.n + net.n*(net.n+net.m-1)).*rand(1 + net.n + net.n*(net.n+net.m-1),1 + net.n + net.n*(net.n+net.m-1));

net.C = [];  %Jacobian Matrix
net.alpha = []; %innovation
net.eta = 0.1;
net.E_not=0;
net.mm = 0;
net.UKF = []; %UKF parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% C++ version has different attributes %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%0.2091*ones(net.m+net.n-1, net.m+net.n-1);
net.c_W = 0.5*(2*rand(net.m+net.n-1, net.m+net.n-1)-1);   %randomly generate values for all possible weights(16x16 in this case)

net.c_W(1:net.m-1,:)=0;                 %zero out 1st 11 rows(no weights between input nodes)       

for i=net.m:net.m+net.n-1               %zero out remaining unused weights
    net.c_W(i,i:end)=0;                 %  taking into account the 'simultaneous' 
end                                     %  weights between the recurrent nodes

%0.00678*ones(net.n,1);
net.c_ww = 0.01*(2*rand(net.n,1)-1);    %randomly generate bias weights

net.c_x = zeros((net.m+net.n-1)*objMaze.size_x*objMaze.size_y,1);
net.store_y = zeros((net.m+net.n-1)*objMaze.size_x*objMaze.size_y,net.MAX_STEPS  );
net.store_f_y = zeros((net.m+net.n-1)*objMaze.size_x*objMaze.size_y,net.MAX_STEPS  );

%% for now make c_W and W identical - only for testing
% % net.c_W = zeros(net.m+net.n-1, net.m+net.n-1);
% % for i=1:net.n
% %     net.c_W(net.m+net.n-i,:) = net.W{i}(net.m+net.n+1-i,2:net.n+net.m);    
% %     net.c_ww(net.n+1-i) = net.W{i}(net.m+net.n+1-i,1);    
% % end

%%% for ALR
net.Lr_W  = 10;
net.Lr_ww = 10;
net.Lr_Ws = 10;
s =(net.n+net.m-1);
net.SUM_F_W_O = ones(s,s);
net.SUM_F_ww_O = ones(net.n,1);
net.SUM_F_Ws_O = 1;

