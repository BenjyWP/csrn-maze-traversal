function [E, c_J, c_DeltaJ]= mazev3_read_output(maze, net, p)
cell_num=0;
for row= 1:maze.size_y
    for col=1:maze.size_x
        cell_num=cell_num+1;
        c_J(row, col) = net.Ws*net.store_y(cell_num*(net.Cell_Nodes-1),p);
        % net.Ws=40; net.store_y= zeros(16*7*7(=784), 10)
        % p=10; net.Cell_Nodes=17;
        % After Feedfoward Calculation, net.store_y stores the output
        % value. Then we take every 16th value, which is the cost function
        % or the ultimate output.
                
        c_DeltaJ(row,col) = c_J(row, col) - maze.A(row, col);
        % c_DeltaJ is the difference between the solution (output) and the 
        % orinal maze.
    end
end
E = sum(sum(c_DeltaJ.*c_DeltaJ));


