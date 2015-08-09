function h1 = drawbox(time, dur, ystart, yheight, color)
% h1 = drawbox(xstart, xlen, ystart, ylen, color);
% 
x = [0 1 1 0]; x = x * dur + time;
y = [0 0 1 1]; y = y * yheight + ystart;

h1 = fill(x,y,color,'FaceAlpha',.5);

return
