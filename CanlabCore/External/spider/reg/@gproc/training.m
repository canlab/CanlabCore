function [dat,algo] =  training(algo,dat)


 disp(['training ' get_name(algo) '.... '])
 [algo.input algo.target]= get_xy(dat);
 [n D] = size(algo.input);
 
 algo.H = zeros(D+2,1);
 [algo.H, fX] = minimize(algo.H, 'gpS00', algo.length, {algo.input, algo.target});
 
 dat=set_x(dat,get_y(dat));
 
 
