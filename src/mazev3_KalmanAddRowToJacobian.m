function net = mazev3_KalmanAddRowToJacobian(net, objMaze, c_DeltaJ,F_W, F_ww, F_Ws)
%%% C++ version

% %% in this version, we add 49 rows to the jacobian, since we consider all
% %% inputs separately
% % OLD STUFF
% % function net = mazev3_KalmanAddRowToJacobian(net, objMaze, DeltaJ,
% c_DeltaJ,F_W, F_ww, F_Ws)
% % net.SUM_F_W=zeros(objMaze.size_y*objMaze.size_x*net.Cell_Nodes,objMaze.size_y*objMaze.size_x*net.Cell_Nodes);   
% % for nn=1:net.MAX_STEPS
% %     net.SUM_F_W = net.SUM_F_W + net.store_f_W{nn};   
% % end
% % net.SUM_F_Ws = net.f_Ws;
% % % now we have 49 outputs, so our Jacobian C is 49 by net.Cell_Nodes^2
% % k = 1;  
% % for  row=1:objMaze.size_y
% %     for col = 1:objMaze.size_x 
% %         % this selects one cell
% %         C_row = [reshape(net.SUM_F_W((k-1)*net.Cell_Nodes+1:k*net.Cell_Nodes,(k-1)*net.Cell_Nodes+1:k*net.Cell_Nodes), 1, net.Cell_Nodes*net.Cell_Nodes)  net.SUM_F_Ws(row, col)];
% %         net.C = [net.C; C_row];               
% %         net.alpha = [net.alpha; DeltaJ(row, col)];
% %         k = k + 1;
% %     end
% % end
% % OLD STUFF ENDS

%----- Total Feedback for c_W weights
s =(net.n+net.m-1);                                     % 5 + 12 - 1 = 16; 
c = objMaze.size_y*objMaze.size_x*s;                    % 7*7*16=784
T_F_W= zeros(net.n,s*objMaze.size_y*objMaze.size_x);    % 5 by 784; 5 recurrent node by all nodes

%F_W = [16 x 7840][16 x (16x49x10 = 7840)]
% c_W = weight matrix = [16x16],  (weights connecting 16 nodes x 16 nodes)...not all are used
% F_W is a feedback matrix. (derivative of output wrt each weight)
% Each cell has its own network therefore
% F_W consist of 49 sets of 16x16 feedback values for these weights (16 x 784)
% since there are 10 recurrent iterations, we have 10 sets of (16 x 784) or [16 x 7840]
% T_F_W takes the sum of the F_W values across the 10 recurrent iterations.
%   or along the rows...across the columns.
% ie. we have 10 sets of 784...so we sum the 1st element in each set
%                                        the 2nd element in each set
%                                        for all 10 sets
% Only rows 12-16 are used...therefore T_F_W only sums along these rows.
for idx=1:s*objMaze.size_y*objMaze.size_x
    %                         12:16 (5 rows)
    %                                          sum each of 784 F_W values(1/col) 
    %                                          across each of the 10 iterations
    T_F_W(:,idx) = sum(F_W(net.m:net.m+net.n-1,idx:c:end),2);    
end

%----- Total Feedback for c_ww weights
c = objMaze.size_y*objMaze.size_x*net.n;                %5x7x7 = 5x49 = 245
T_F_ww = zeros(1,net.n*objMaze.size_y*objMaze.size_x);  %T_F_ww = [1 x 245] ( 1 x 5*49)

%F_ww = 1 by 2450.....10 sets of 1x245
%sum each of 245 F_ww values(1/col) across each of the 10 iterations
for idx=1:net.n*objMaze.size_y*objMaze.size_x 
    T_F_ww(:,idx) = sum(F_ww(idx:c:end));    
end

k=1;
for  row=1:objMaze.size_y
    for col = 1:objMaze.size_x 
      % this selects one cell
      % s=16; k is varing, T_F_W - 5 by 784 but now reshaped into 1 by 86
        C_row = [reshape(T_F_W(:,(k-1)*s+1:k*s), 1, net.n*s)  T_F_ww(:,(k-1)*net.n+1:k*net.n) F_Ws(k)];
                % 5 by 16 into 1 by 80    +       1 by 5 +  1 = 86 (F_Ws: 1 by 49)         
        net.C = [net.C; C_row];               
        net.alpha = [net.alpha; c_DeltaJ(row, col)]; %(net.alpha: 1x49)
        k = k + 1;
    end
end

size(T_F_W);