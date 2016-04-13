function mov = make3Davi_uncompressed(varargin)
% :Usage:
% ::
%
%    mov = make3davi([opt] Options_Structure)
%
% Makes an avi movie file called head3d[x].avi or whatever you specify.
%
% :Options:
%
%   **O.name:**
%        'output_name.avi';
%
%   **O.fps:**
%        frames per second
%
%   **O.length:**
%        length of movie in s
%
%   **O.H:**
%        cluster handles in vector   (for adjusting transparency with time)
%
%   **O.timecourse{i}:**
%        cell array of time courses for each cluster
%
%   **O.timeres:**
%        resolution, in s, of timecourse data
%
%   **O.azOffset:**
%        azimuth value to offset, positive = move clockwise
%
%   **O.elOffset:**
%        elevation value to move through, positive = inf to superior
%
%   **O.zoom:**
%        zoom value to end up with
%
%   **O.add2movie:**
%        add to existing movie - enter mov structure in this field
%
%   **O.closemovie:**
%        1 or 0, close the movie afterward or not.
%
% Blank fields for az, el, zoom, timecourse indicate that these functions should not be performed
% This script spirals up, right, and in 36 degrees
%
% :Notes: my indeo5 one wouldn't work in media player
%         also only seems to work if you add clusters before head isosurfaces
%
% Start with the image in the location you want to zoom in on
% but with no zoom.
%
% ..
%    By Tor Wager, 10/07/01
% ..


% ..
%    set default values and inputs
% ..
myName = 'head3d1.avi'; i = 1;                      % set name to first unused file #
while exist(myName) == 2, myName = ['head3d' num2str(i) '.avi'];, i = i+1;,end
O.name = myName;
O.fps = 10;
O.length = 8;
O.H = [];
O.timecourse = [];
O.timeres = 1;
O.azOffset = 36;                                    % positive starts from left, moves r
                                                    % spins brain clockwise
O.elOffset = 20;                                    % positive starts from below, moves up
O.zoom = 2;                                         % final image will be this many times closer than at start
O.add2movie = [];                                   % existing movie structure to add to.  [] = new movie file.
O.closemovie = 1;                                   % close movie after finishing.  0 to add more frames later.

if nargin > 0
    inO = varargin{1};
    if isfield(inO,'name'), O.name = inO.name;,end
    if isfield(inO,'fps'), O.fps = inO.fps;,end
    if isfield(inO,'length'), O.length = inO.length;,end
    if isfield(inO,'H'), O.H = inO.H;,end
    if isfield(inO,'timecourse'), O.timecourse = inO.timecourse;,end
    if isfield(inO,'timeres'), O.timeres = inO.timeres;,end
    if isfield(inO,'azOffset'), O.azOffset = inO.azOffset;,end
    if isfield(inO,'elOffset'), O.elOffset = inO.elOffset;,end
    if isfield(inO,'zoom'), O.zoom = inO.zoom;,end
    if isfield(inO,'add2movie'), O.add2movie = inO.add2movie;,end
    if isfield(inO,'closemovie'), O.closemovie = inO.closemovie;,end
end



numFrames = O.fps .* O.length;
interpValue = O.fps ./ O.timeres;

% -------------------------------------------
% * check values to make sure everything's ok
% -------------------------------------------

if ~isempty(O.timecourse) & length(O.H) ~= length(O.timecourse), error('Number of handles and timecourses must be equal.'),end
if interpValue ~= round(interpValue), error('Frames per second must be a multiple of time res - adjust either one.'),end
if any(~ishandle(O.H)), 
    warning('At least one cluster handle is not a valid object. Searching for handles by Tag name cluster[x]')
    for i = 1:length(O.H)
        O.H(i) = findobj('Tag',['cluster' num2str(i)]);
    end
    if any(~ishandle(O.H)), error('Cannot find valid cluster handles.'),end
end

% -------------------------------------------
% * create transparency values from timecourses
% -------------------------------------------
% normalize all timecourses to alpha values, from 0 to 1, by dividing by max value
if ~isempty(O.timecourse)
  for i = 1:length(O.timecourse), 
    %O.timecourse{i} = O.timecourse{i} - min(O.timecourse{i});       % prevent negative values, set min to 0
    %min(O.timecourse{i})
    %O.timecourse{i} = O.timecourse{i} ./ max(abs(O.timecourse{i})); % max is one
    numSeconds = length(O.timecourse{i}) .* O.timeres;
    alphaVal{i} = interp(O.timecourse{i},interpValue);              % interpolate to number of frames
    alphaVal{i} = alphaVal{i} - min(alphaVal{i});                   % prevent negative values, set min to 0
    alphaVal{i} = alphaVal{i} ./ max(abs(alphaVal{i}));             % max is one
    alphaVal{i} = [ alphaVal{i} zeros(1,min([0 numFrames - length(alphaVal{i})])) ];   % pad with zeros
    
    if any(alphaVal{i} < 0 | alphaVal{i} > 1),
        alphaVal{i},min(alphaVal{i}),max(alphaVal{i})
        error('Script error: alpha values must be between 0 and 1.')
    end
  end
end



% -------------------------------------------
% * set initial view and parameters
% -------------------------------------------
set(gcf,'Color','k')
axis off
[az,el] = view;
if ~isempty(O.azOffset),az = az - O.azOffset;,az2add = O.azOffset ./ numFrames;,end
if ~isempty(O.elOffset),el = el - O.elOffset;,el2add = O.elOffset ./ numFrames;,end
view(az,el)
lightFollowView
if ~isempty(O.zoom), myZoom = 1 + ((O.zoom-1) ./ numFrames);, end
   


% -------------------------------------------
% * initialize movie file, if necessary
% -------------------------------------------

if isempty(O.add2movie)
    try
        mov = avifile(O.name,'Quality',100,'Compression','None','Fps',O.fps);
    catch
        error('Cannot open movie file.  Probably already open: Try changing filename.')
    end
else 
    mov = O.add2movie;
end

%try
    
H = gca;

% -------------------------------------------
% * rotate, zoom, etc. and add frames
% -------------------------------------------
for i = 1:numFrames
    
	if ~isempty(O.azOffset),az = az + az2add;,end
    if ~isempty(O.elOffset),el = el + el2add;,end
	if ~isempty(O.azOffset) | ~isempty(O.elOffset)
        view(az,el);
        lightFollowView
    end
	if ~isempty(O.zoom),camzoom(myZoom);,end
    if ~isempty(O.timecourse)
        for j = 1:length(alphaVal)
            set(O.H(j),'FaceAlpha',alphaVal{j}(i))
        end
        myTime = i ./ O.fps;
        if myTime == round(myTime)
            title(['Time = ' num2str(i ./ O.fps) '.0 s'],'Color','w','FontSize',18)
        else
            title(['Time = ' num2str(i ./ O.fps) ' s'],'Color','w','FontSize',18)
        end
        
    end
    drawnow
	try
        mov = addframe(mov,H);
    catch
        error('Cannot write frame.  Failed to set stream format??')
    end
    
end


%catch, warning('There was an error making this movie file.  No frames written.')
%end

% -------------------------------------------
% * close the movie
% -------------------------------------------

if O.closemovie, mov = close(mov);, end


return
