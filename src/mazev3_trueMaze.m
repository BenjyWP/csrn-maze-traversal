function [maze] = mazev3_trueMaze(input_maze)
maze = input_maze;
[r,c] = find(maze==1);
maze(r,c) = 0;
s = 0;
while sum(sum(maze)) ~= s
    s = sum(sum(maze));    
    for r=2:size(input_maze,1)-1
        for c=2:size(input_maze,2)-1
            if and(input_maze(r,c)~=25 , input_maze(r,c)~=1)
                maze(r,c) = 1 + min(maze(r-1,c), min( min(maze(r+1,c),maze(r,c-1)), maze(r,c+1)));
            end
        end
    end
end