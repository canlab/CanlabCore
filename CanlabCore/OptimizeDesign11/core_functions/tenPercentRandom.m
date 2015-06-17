function lM = tenPercentRandom(lM,saveWhich,conditions)
% function lM = tenPercentRandom(lM,saveWhich,conditions)
%
% saves the best list from randomization, and selects one
% of the condition numbers at random for 10% of the list 
% values across the population.

    %conditions = [conditions 0];
    mybest = lM(:,saveWhich);
    myfirst = lM(:,1);
    mylen = length(conditions);
    for i = 1:.1 * size(lM,1) .* size(lM,2)
        a = ceil(rand*mylen);
        %a = randperm(mylen);
        %a = a(1);
        a = conditions(a);
        
        b = ceil(rand*size(lM,1));
        c = ceil(rand*size(lM,2));
        %b = randperm(size(lM,1)); b = b(1);
        %c = randperm(size(lM,2)); c = c(1);
        lM(b,c) = a;
    end
    lM(:,saveWhich)  = mybest;
    lM(:,1) = myfirst;
return