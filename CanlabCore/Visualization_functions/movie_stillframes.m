function mov = movie_stillframes(numframes,mov) 
% Add still frames to a movie

if nargin < 2 || isempty(mov)
    mov = avifile('mymovie.avi','Quality',75,'Compression','None','Fps',8);
end

 axis vis3d

 for i = 1:numframes
     
    H = gca;
    drawnow

    try
    mov = addframe(mov,H);
    catch
    disp('Cannot write frame.  Failed to set stream format??')
    mov = close(mov);
    end
 end
 
    
return
