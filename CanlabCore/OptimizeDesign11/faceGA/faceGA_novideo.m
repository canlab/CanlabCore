function bestp = faceGA
%
% A silly face example of a genetic algorithm.
% By Tor Wager, June 2002
%
% calls functions:
% drawFace.m
% breedorgs.m   (from OptimizeDesign10)

gensize = 100;
numgens = 100;
numopt = 17;
change_every = 30;

% ------------------------------------------
% * set up figures
% ------------------------------------------

fgen = figure;set(gcf,'MenuBar','none');set(gcf,'Name','Faces in This Generation')
set(gcf,'Position',[30          70         975        1040])
set(gcf,'NumberTitle','off');set(gcf,'Color','w')

% for movie single figure
for i = 1:30, h(i) = subplot(6,5,i);, mypos(i,:) = get(h(i),'Position');,end
mypos(:,1) = mypos(:,1) - .1;
mypos(:,3) = mypos(:,3) * .7;
mypos(:,1) = mypos(:,1) * .7;
clf
for i = 1:30, nh(i) = axes('position',mypos(i,:));, end
set(gcf,'Position',[30          70        1514        1040])
hgen = axes('position',[.591 .5347 .36 .36]);
hideal = axes('position',[.591 .11 .36 .36]);

%mov = avifile('facega.avi','Quality',75,'Compression','none','Fps',10);

%bgen = figure;set(gcf,'MenuBar','none');set(gcf,'Name','Best of This Generation')
%set(gcf,'Position',[1021         719         516         393])
%set(gcf,'NumberTitle','off');set(gcf,'Color','w')
%tgen = figure;set(gcf,'MenuBar','none');set(gcf,'Name','Target Face')
%set(gcf,'Position',[1024         195         516         393])
%set(gcf,'NumberTitle','off');set(gcf,'Color','w')

% ------------------------------------------
% * set up lists
% ------------------------------------------
p = rand(gensize,20);
idealface = [1.0000    1.0000    1.0000    1.0000         0    0.4486    0.5244    0.1715    0.5000    1.0000    0.2000    0.1414    0.4570    1.0000 0.8000    0.2248    0.9089    0.0073    0.5887    0.5421]; 
idealfit = repmat(idealface,gensize,1);
%figure(tgen); drawFace(idealface)
axes(hideal); cla; drawFace(idealface); title('Face To Match','FontSize',14)

for gen = 1:numgens
    
    if mod(gen,change_every) == 0
        idealface = rand(1,20);
        idealfit = repmat(idealface,gensize,1);
        axes(hideal); cla; drawFace(idealface)
    end
    
    for g = 1:gensize
        
        if g < 31 
            
            axes(nh(g))
            %figure(fgen)
            %subplot(6,5,g)
            drawFace(p(g,:))
            %H = gcf; mov = addframe(mov,H);
        end
        
    end
    
    % compute fitness and save best
    f = -1 * (sum(((idealfit(:,1:numopt) - p(:,1:numopt)) .^ 2),2));
    bestp = p(f == max(f),:); bestp = bestp(1,:);
    %figure(bgen); clf; drawFace(bestp)
    axes(hgen); cla; drawFace(bestp); title('Best of This Generation','FontSize',14)
    
    fm = f(f == max(f));
    text(.8,.8,1,num2str(round(fm(1)*1000)),'FontSize',18)
    
    f = (f - mean(f)) ./ std(f);
    
    % crossbreed
    p = breedorgs(f',p',2.1)';

    H = gcf; %mov = addframe(mov,H);
end
        
%mov = close(mov);  