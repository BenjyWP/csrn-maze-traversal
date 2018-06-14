function [pct] = mazev3_goodness_of_solution(maze, solution)
% function [pct] = mazev3_goodness_of_solution(maze, solution)
% returns what percentage of the maze cells has the right direction
% maze is the true solution
% solution is the approximate solution being tested
%% agree on directions
%% 0 - STAY, means no where to go
%% 1 - UP, 2 - LEFT, 3-DOWN, 4 - RIGHT
OBST = max(max(maze));
pct = 0;
APPORX_DIR = getDirmatrix(maze, solution);%The directions of the solutions
EXACT_DIR = getDirmatrix(maze, maze); %The directions of the maze
%% determine what pct is correct
correct = max(size(find(APPORX_DIR(2:end-1,2:end-1) -EXACT_DIR(2:end-1,2:end-1) ==0)));
cells_not_obst = max(size(find(maze(2:end-1,2:end-1)~=OBST )));
cells_obst = max(size(find(maze(2:end-1,2:end-1)==OBST )));
pct = 100*(correct-cells_obst)/cells_not_obst;

function dirm = getDirmatrix(maze, m)
OBST = max(max(maze));
dirm = zeros(size(m));
mindir = 0;
minval = 0;
for row = 2:size(m, 1)-1
    for col = 2:size(m, 2)-1
        if maze(row,col) ~= OBST
            mindir = 0;
            minval = 0;
            if m(row, col)>m(row-1,col)
                mindir = 1;
                minval = m(row-1,col);
            else
               minval = m(row,col); 
            end
            if and(m(row, col)>m(row,col-1), m(row,col-1)<minval)
                mindir = 2;
                minval = m(row,col-1);
            end
            if and(m(row, col)>m(row+1,col), m(row+1,col)<minval)
                mindir = 3;
                minval = m(row+1,col);
            end
            if and(m(row, col)>m(row,col+1), m(row,col+1)<minval)
                mindir = 4;
                minval = m(row,col+1);
            end 
            dirm(row, col) = mindir;
        end
    end
end
