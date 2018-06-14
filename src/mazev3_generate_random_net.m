function [best_net, MIN_EE] = mazev3_generate_random_net(MAZES, TRIALS )
MIN_EE=Inf;
['Generating Random Nets...']
for k=1:TRIALS
    display([num2str(k) ' ... ']);
    maze=MAZES{1};
    net =  mazev3_initialize_net(maze);
    EE=0;
    for maze_num=1:size(MAZES,2)
        maze=MAZES{maze_num};
        p=net.MAX_STEPS;
        net = mazev3_assign_external_inputs(maze, net);
        
        [y, v, net.store_y]=MAZE_CELLULAR_SRN(net.c_W, net.c_ww, net.c_x, net.n, net.m-1, maze.size_x, p,0,0,0,0);
        [E , c_J, c_DeltaJ]= mazev3_read_output(maze, net, p);
        EE = EE+E;    
    end  %MAZES
    
    if EE < MIN_EE
        MIN_EE = EE;
        best_net = net;
    end
end  % trials    
['Done. Best Error is ' num2str(MIN_EE)]