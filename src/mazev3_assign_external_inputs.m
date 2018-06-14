function net = mazev3_assign_external_inputs(maze, net)
% assign initial values for recurent links and target/goal inputs here
net.c_x = zeros(net.m+net.n-1,1);                       %create input vector for 1st cell
net.c_x(3:3+4+net.n-2)=-1;                              %initial neighbor and self recurrent inputs
net.c_x = repmat(net.c_x, maze.size_x*maze.size_y,1);   % replicate for other cells
                                                        % net.c_x  49 by 1 times
                                                        % net.c_x =16; 16*49 =784

% use the target map to set the values for the 1st 2 inputs(x1 & x2) for each cell 
for row=1:maze.size_y
    for col=1:maze.size_x
        c = (row-1)*maze.size_x + col;                  %every cell : 1-49
        if maze.B(row, col) == 25
            net.c_x((c-1)*(net.Cell_Nodes-1) + 1) = 1;  %maze.B(row, col) = target value
        end                                             %every 16th value is re-evaluated 
        if maze.B(row, col) <= 1
	        net.c_x((c-1)*(net.Cell_Nodes-1) + 2) = 1;  %target goal values
        end                                             %the only 17th value that has 1 
    end
end                                                     % eventually we have 36 of value '1' 