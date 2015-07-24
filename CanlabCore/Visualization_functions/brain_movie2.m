
% set up brain

p = addbrain('brainbottom');

color = get(p(1),'FaceColor');
set(p(1),'Visible','off');
set(p(2),'Visible','off');

ps = addbrain('transparent_surface');
set(ps,'FaceColor',[.6 .4 .3]);
set(ps,'FaceAlpha',1);



pa = addbrain('amygdala');
ph = addbrain('hippocampus');

pm = addbrain('midbrain');
pt = addbrain('thalamus');
pn = addbrain('nucleus accumbens');

pc = addbrain('caudate');
pp = addbrain('putamen');
pg = addbrain('globus pallidus');

pn = addbrain('nucleus accumbens');

axis image; axis tight; axis vis3d

material dull
camzoom(1.3);



mov = avifile('mymovie2.avi','Quality',75,'Compression','None','Fps',5);

% opaque surface
view(225,0)
set(ps,'FaceAlpha',1);

[az,el]=view;

myaz = linspace(az,135,20);
myel = linspace(el,30,20);
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
    lightFollowView,

end

% add activations


% make transparent

set(p(1),'Visible','on');
set(p(2),'Visible','on');
for i = linspace(1,.05,12)

    set(ps,'FaceAlpha',i);

    drawnow
    try
        mov = addframe(mov,H);
    catch
        disp('Cannot write frame.  Failed to set stream format??')
        mov = close(mov);
    end

end

[az,el]=view;

myaz = linspace(az,225,20);
myel = linspace(el,20,20);
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
    lightFollowView,

end






