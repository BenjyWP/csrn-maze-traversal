function [maz, maze_j] = mazev3_generate_random_maze(maze_size, max_obst)
%function [maze, maze_j] = mazev2_generate_random_maze(maze_size, max_obst)
%returns 7 by 7 matrices corresponding to A and B
%what if the generated maze contains unreachable cells - 
% their value will be obstacle+1  
obstacle = 25;
goal = 1; 
maz = ones(maze_size+2,maze_size+2)*24;
maz(1,:)=obstacle;maz(end,:)=obstacle;
maz(:,1)=obstacle;maz(:,end)=obstacle;
taken_pos = [];
r = floor(rand()*maze_size+ 1);
c = floor(rand()*maze_size+ 1);
maz(r+1,c+1) = goal;
taken_pos = [taken_pos; r c];
obsts = maze_size + floor(rand()*(max_obst-maze_size)+ 1);
for i=1:obsts
    f = [1];
    while ~isempty(f) 
        r = floor(rand()*maze_size+ 1);
        c = floor(rand()*maze_size+ 1);    
        f = find(and(taken_pos(:,1)==r, taken_pos(:,2)==c));
        maz(r+1,c+1) = obstacle;
        taken_pos = [taken_pos; r c];
    end
end

maze_j = mazev3_trueMaze(maz);
%maze_v2_display_maze(maze_j,1);