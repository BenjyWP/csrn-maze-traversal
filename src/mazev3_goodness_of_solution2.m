%% goodness of navigation
%% count from how many cells the goal can be reached
X = MAZES{1}.A;
for r = 1:size(X, 1)
    for c=1:size(X,2)
        if X(r,c)~=25
            currPosRow = r; currPosCol=c; Cost = 0;
            while (X(r,c)~=0)
                % up, left, down, right
                [v i]=min([X(r-1,c) X(r,c-1) X(r+1,c) X(r,c+1) X(r,c)]);
                switch i
                    case 1
                        currPosRow=currPosRow-1; Cost=Cost+1;
                    case 2
                        currPosCol=currPosCol-1; Cost=Cost+1;
                    case 3
                        currPosRow=currPosRow+1; Cost=Cost+1;
                    case 4
                        currPosCol=currPosCol+1; Cost=Cost+1;
                    case 5
                        if X(currPosRow,currPosCol)==0
                            %goal reached
                        else
                            %local sink, dead end
                        end
            end
        end
    end
end