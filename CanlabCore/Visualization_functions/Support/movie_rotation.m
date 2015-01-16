function mov = movie_rotation(targetaz,targetel,mov) 
% mov = movie_rotation(targetaz,targetel,mov)
%
% Make a movie of a rotating brain
%
% Examples:
%
% mov = movie_rotation(135,30)
% mov = movie_rotation(250,30,mov)
% 
% Batch mode:
% mov = movie_rotation('batch360.1');
% mov = movie_rotation('right2left');
%
% tor wager, may 06

% Batch-mode: pre-set movies
if strcmp(targetaz,'batch360.1')
    axis vis3d
    axis image
    view(0,90); axis off; set(gcf,'Color','w');
    mov = movie_rotation(135,30);
    mov = movie_rotation(250,30,mov);
    mov = movie_rotation(360,90,mov);
    mov = close(mov);
    return
elseif strcmp(targetaz,'right2left')




if nargin < 3 || isempty(mov)
    mov = avifile('mymovie.avi','Quality',90,'Compression','None','Fps',10);
end

 % add to existing
 %O = struct('add2movie',[],'zoom',1,'azOffset',[],'elOffset',[],'timecourse',[],'fps',5,'length',6);
 
 axis vis3d
 
 [az,el]=view;
 
 myaz = linspace(az,targetaz,20);
 myel = linspace(el,targetel,20);

 for i = 1:length(myaz)
     
    H = gca;
    drawnow

    try
    mov = addframe(mov,H);
    catch
    disp('Cannot write frame.  Failed to set stream format??')
    mov = close(mov);
    end

    view(myaz(i),myel(i));

    axis image
    lightRestoreSingle(gca);
    %     if i >= 10, lightFollowView, end
    %     
    %     if i == 10, 
    %         lighting phong
    %     end
    
 end
    
     try
    mov = addframe(mov,H);
    catch
    disp('Cannot write frame.  Failed to set stream format??')
    mov = close(mov);
     end
    
return
    
 