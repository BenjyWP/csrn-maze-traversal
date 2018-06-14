function [ ret, override] = mazev3_stoppingcondition(net,maze, ERRORS)
k=1;
override = 0;
count = 0;
if size(ERRORS,1)<2*k+1
    ret = 0;
    return;
end
if ERRORS(end,3) > ERRORS(end-k,3) 
    count = count+ 1;
end
if ERRORS(end-k,3) > ERRORS(end-2*k,3) 
    count = count + 1;
end
if count == 2
    ret = 1;
else
    ret = 0;
end

if ret==1
    reply = input('Stop Now (V - never)? Y/N/V [Y]: ','s');
    if isempty(reply)
        reply = 'Y';
    end
    if reply=='N'
        ret=0;
    end
    if reply=='V'
        ret=0;
        override = 1;
    end
end