function [ output_args ] = maze_v3_display_maze( J , new_figure)
if new_figure == 1
    figure;    
    clf;
    axis;
end
for x=1:size(J,1)-2
    for y=1:size(J,2)-2
        rectangle('Position', [x,y,1,1]);
    end
end
for x=1:size(J,1)-2
    for y=1:size(J,2)-2
        text(y+0.3, size(J,1)-1-x+0.5, num2str(J(x+1,y+1)))
    end
end