fn = fopen(filename,'a');
fprintf(fn, '\n');
fprintf(fn, '%d;', max(size(MAZES)));
fprintf(fn, '%d;', experiment_counter);
fprintf(fn, '%d;', MAZES{1}.size_x-2);
fprintf(fn, '%d;', t);
fprintf(fn, '%d;', EE);
fprintf(fn, '%d;', NUM_RAND_NETS);
fprintf(fn, '%s;', net.METHOD);
fprintf(fn, '%s;', datestr(now));
fprintf(fn, '%s;', num2str(net.n));
fprintf(fn, '%s;', num2str(net.N));
fprintf(fn, '%s;', formula(net.inline_r));
fprintf(fn, '%d;', TOT_E);
fprintf(fn, '%d;', best_t);
fprintf(fn, '%d;', TOT_PCT/size(MAZES,2));
fprintf(fn, '%d;', TOT_TEST_PCT/size(MAZES,2));
fclose(fn);


res = dir('maze_struct.mat');
if size(res,1)==0
    maze_struct = cell(1,3);
else    
    load maze_struct
end
cc = size(maze_struct,1)+1;
maze_struct{cc,1} = experiment_counter;
maze_struct{cc,2} = MAZES;
maze_struct{cc,3} = TEST_MAZES;
save maze_struct.mat maze_struct
