function [MAZES] = mazev3_generate_random_mazes(num_mazes, maze_size, max_obst, STR_TITLE)
MAZES = cell(1,1);
num = 0;
while num  < num_mazes 
    [maz, maze_j] = mazev3_generate_random_maze(maze_size, max_obst);
    while max(max(maze_j)) > 25
        [maz, maze_j] = mazev3_generate_random_maze(maze_size, max_obst);
    end
    num = num + 1;
    objMaze.size_x =maze_size+2;
    objMaze.size_y =maze_size+2;
    objMaze.B = maz;
    objMaze.A = maze_j;
    MAZES{num} = objMaze;    
end
return
figure
if num_mazes <= 4
    for i=1:num_mazes
        subplot(2,2,i)
        if i==1 
            title(STR_TITLE);
        end
        maze_v3_display_maze(MAZES{i}.A,0);                
        axis([1 maze_size+1 1 maze_size+1]);
    end
elseif num_mazes <= 6
    for i=1:num_mazes
        subplot(3,2,i)
        if i==1 
            title(STR_TITLE);
        end
        maze_v3_display_maze(MAZES{i}.A,0);                
        axis([1 maze_size+1 1 maze_size+1]);
    end
elseif num_mazes <= 9
    for i=1:num_mazes
        subplot(3,3,i)
        if i==1 
            title(STR_TITLE);
        end
        maze_v3_display_maze(MAZES{i}.A,0);
        axis([1 maze_size+1 1 maze_size+1]);
    end
elseif num_mazes <= 12
    for i=1:num_mazes
        subplot(3,4,i)
        if i==1 
            title(STR_TITLE);
        end        
        maze_v3_display_maze(MAZES{i}.A,0);                
        axis([1 maze_size+1 1 maze_size+1]);
    end
elseif num_mazes <= 16
    for i=1:num_mazes
        subplot(4,4,i)
        if i==1 
            title(STR_TITLE);
        end        
        maze_v3_display_maze(MAZES{i}.A,0);                
        axis([1 maze_size+1 1 maze_size+1]);
    end
end

%maze_v2_display_maze(maze_j,1);