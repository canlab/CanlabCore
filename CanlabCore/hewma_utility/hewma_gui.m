function varargout=hewma_gui(Action,varargin)
% :Usage:
% ::
%
%     varargout=hewma_gui(Action,varargin)
%
% To run, type "hewma" at the Matlab prompt
%
% ..
%    Tor Wager 
%
%    Thanks to Tom Nichols for the excellent GUI shell!
% ..


global EXPT % global variables we need for this shell
global cl
    
%-Format arguments
%-----------------------------------------------------------------------
if nargin == 0, Action='Init'; end


switch lower(Action)
    
    case lower('Init')
    %=======================================================================

    %clc
    %BrainVowager_defaults;
    hewma_gui('AsciiWelcome')
    hewma_gui('CreateMenuWin')

    % load DB

    if ~isempty(EXPT)
        disp('Using EXPT already in memory.'); 
    elseif exist('EXPT.mat') == 2
        disp('loading EXPT.mat: Data structure (in current dir) is available.'); load EXPT;  
        
    else
        disp('You need to load or create an EXPT.mat file containing experiment information ');
        disp('for full functionality.  Click Setup or see CreateExpt.m')
        disp('(No EXPT.mat file in current directory.)');
        EXPT = [];
        fprintf(1,'\n')
    end
        
    varargout{1} = EXPT;    
        
    case lower('AsciiWelcome')
    %=======================================================================
    disp( 'Welcome to the hewma change-point analysis gui.  HEWMA was written by Martin Linquist and Tor Wager, Jan 2006')
    fprintf('\n')

    case lower('Ver')
    %=======================================================================
    varargout = {'Hewma Menu'};
  


    case lower('CreateMenuWin')
    %=======================================================================
    close(findobj(get(0,'Children'),'Tag','hewma_gui Menu'))

    
    %-Initialize hewma_gui menu window
    %-----------------------------------------------------------------------
    [F, winwid, winh] = hewma_gui('initFigure');
    
    
    % default button sizes and positions, etc.
    
    topbutton = winh-100;        % y location of top button
    butspace = 30;               % y spacing of buttons
    
    fullbutxy = [160 25];       % full-length button width and height
    halfbutxy = [80 25];        % (left-hand) half-width button w and h
    rightbutxy = halfbutxy;
    rightbutx = 110+halfbutxy(1)+5;  % right-hand button start x
    
    popupxy = [135 25];
    goxy = [20 25];
    



    %-Frames and text
    %-----------------------------------------------------------------------
    axes('Position',[0 0 80/winwid winh/winh],'Visible','Off')
    text(0.5,0.475,'HEWMA Toolbox',...
        'FontName','Times','FontSize',36,...
        'Rotation',90,...
        'VerticalAlignment','middle','HorizontalAlignment','center',...
        'Color',[1 1 1]*.6);

    text(0.2,0.96,'Change point fMRI analysis',...
        'FontName','Times','FontSize',16,'FontAngle','Italic',...
        'FontWeight','Bold',...
        'Color',[1 1 1]*.6);

    uicontrol(F,'Style','Frame','Position',[095 005 winwid-100 winh - 30],...
        'BackgroundColor',hewma_gui('Color'));  % colored frame
    uicontrol(F,'Style','Frame','Position',[105 015 winwid-120 winh - 50]);  % inner gray frame

    %-Buttons to launch hewma_gui functions
    %-----------------------------------------------------------------------
    
    % -------------------------------------------
    % Section - EWMA (single subject)
    uicontrol(F,'Style','Text',...
        'String','EWMA (first-level)','FontSize',14,...
        'HorizontalAlignment','Center',...
        'Position',[115 topbutton+30 fullbutxy],...
        'ForegroundColor','y','FontWeight','b');
    % -------------------------------------------
    
    % -------------------------------------------
        
    % EWMA setup
    buttontext = 'Setup';
    str = 'CreateExpt(''hewma'');';                          % callback function
    str = hewma_gui('ExpandString',str);                 % display then execute
    uicontrol(F,'String',buttontext,...
        'Position',[110 topbutton-(1-1)*butspace halfbutxy],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');
    
    % EWMA analysis
    buttontext = 'Run EWMA';
    str = 'EXPT = wb_multisubject_ewma(EXPT);';             % callback function
    str = hewma_gui('ExpandString',str);                 % display then execute
    uicontrol(F,'String',buttontext,...
        'Position',[rightbutx topbutton-(1-1)*butspace rightbutxy],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');
    
    
        
    % -------------------------------------------
    % Section - HEWMA (multi-subject)
    uicontrol(F,'Style','Text',...
        'String','HEWMA (group)','FontSize',14,...
        'HorizontalAlignment','Center',...
        'Position',[115 topbutton-(2-1)*butspace fullbutxy],...
        'ForegroundColor','y','FontWeight','b');
    % -------------------------------------------
    
    % -------------------------------------------
        
    % HEWMA setup
    buttontext = 'Setup';
    str = 'disp(''Done, if you ran EWMA first and saved EXPT.'');';                          % callback function
    str = hewma_gui('ExpandString',str);                 % display then execute
    uicontrol(F,'String',buttontext,...
        'Position',[110 topbutton-(3-1)*butspace halfbutxy],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');
    
    % HEWMA analysis
    buttontext = 'Run HEWMA';
    str = 'EXPT = wb_hewma_shell(EXPT,1,0,EXPT.mask);';             % callback function
    str = hewma_gui('ExpandString',str);                 % display then execute
    uicontrol(F,'String',buttontext,...
        'Position',[rightbutx topbutton-(3-1)*butspace rightbutxy],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');
    
        
 
    % -------------------------------------------
    % Section - HEWMA Results
    uicontrol(F,'Style','Text',...
        'String','Results','FontSize',14,...
        'HorizontalAlignment','Center',...
        'Position',[115 topbutton-(4-1)*butspace fullbutxy],...
        'ForegroundColor','y','FontWeight','b');
    % -------------------------------------------
    
    % -------------------------------------------
        
    % HEWMA results -- activation
    buttontext = 'Activation';
    str = 'hewma_gui(''activation'',''hewma_sig.img'');';                          % callback function
    %str = hewma_gui('ExpandString',str);                         % display then execute
    uicontrol(F,'String',buttontext,...
        'Position',[110 topbutton-(5-1)*butspace halfbutxy],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');
  
    % HEWMA results -- activation popup menu
    buttontext = {'Map Activation' 'Map Group Diffs' 'Map Onset Time' 'Map Duration'};
    
    str = 'cl = hewma_gui(''activation'');';                          % callback function
    %str = hewma_gui('ExpandString',str);                         % display then execute
    pop1 = uicontrol(F,'Style','popupmenu', 'Tag','ImagePop', ...
    'String',buttontext,...
        'Position',[110 topbutton-(5-1)*butspace fullbutxy],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');
 
    % HEWMA results -- timeseries checkbox
    buttontext = {'Timeseries plot'};

    uicontrol(F,'Style','checkbox', 'Tag','TimeseriesChk', ...
    'String',buttontext,...
        'Position',[110 topbutton-(6-1)*butspace fullbutxy],...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');
            
               
    % HEWMA results -- classification popup menu
    buttontext = {'Classify CP' 'Classify CP + Dur'};
    
    str = 'cl = hewma_gui(''classify'');';                          % callback function
    %str = hewma_gui('ExpandString',str);                         % display then execute
    pop1 = uicontrol(F,'Style','popupmenu', 'Tag','ClassifyPop', ...
    'String',buttontext,...
        'Position',[110 topbutton-(7-1)*butspace fullbutxy],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');
        

    % Blob Display Tool 
    if exist('scn_roi_gui') == 2
        str = 'scn_roi_gui;';                                           % callback function
    else
        str1 = 'if isempty(cl), disp([''Load a clusters cl.mat file first!'']), return, end;';
        str2 = 'cluster_orthviews(cl,{[1 0 0]}); set(gcf,''WindowButtonUpFcn'',''Meta_interactive_table;'')';
        str = [str1 str2];
    end
    
    str = hewma_gui('ExpandString', str);                   % display then execute
    uicontrol(F,'String','Blob Display Tool',...
        'Position',[110 topbutton-(8-1)*butspace fullbutxy],...
        'CallBack',str,...
        'Interruptible','on',...
        'ForegroundColor','k','FontWeight','b');    
   
    
    set(F,'Pointer','Arrow','Visible','on')



    case lower('Color')
    %=======================================================================
    % hewma_gui('Color')
    %-----------------------------------------------------------------------
    % %-Developmental livery
    % varargout = {[0.7,1.0,0.7], 'Lime Green'};
    %-Distribution livery
    varargout = {[.2 0.2 .6], 'Purple'};

    
    case lower('ExpandString')
    %=======================================================================
    % hewma_gui('ExpandString')
    % Expand an action button callback string (a command string to be
    % evaluated)
    % so that it first displays the command, and then executes it
    %-----------------------------------------------------------------------      
        str = varargin{1}; str2 = [];
        for i = 1:length(str)
            if str(i) == '''', str2(end+1) = ''''; str2(end+1) = '''';
            else, str2(end+1) = str(i);
            end
        end

        str = ['disp(''' char(str2) '''), ' str ];     % display then execute

        varargout = {str};
  
        
    case lower('initFigure')
    %=======================================================================
    % [F, winwid, winh] = hewma_gui('initFigure')
    %-----------------------------------------------------------------------
    % Get the position of the main BrainVowager menu, or if
    % not available, default screen pos.
        
        % default sizes, etc.
    S = get(0,'ScreenSize');

    winwid = 300;               % window width
    winh = 400;                 % window height
    pos = [S(3)/2+150,S(4)/2-140,winwid,winh];  % default
    
    h = findobj('Tag','hewma_gui Menu');
    if ~isempty(h), 
        pos = get(h,'Position'); 
        winwid = pos(3); winh = pos(4);
        pos(1) = pos(1) + winwid;   % put next to main figure
    end
    
    %-Open hewma_gui menu window
    %----------------------------------------------------------------------
    
    F = figure('Color',[1 1 1]*.8,...
        'Name',hewma_gui('Ver'),...
        'NumberTitle','off',...
        'Position',pos,...
        'Resize','off',...
        'Tag','hewma_gui Menu',...
        'Pointer','Watch',...
        'MenuBar','none',...
        'Visible','off');
    
    varargout{1} = F; varargout{2} = winwid; varargout{3} = winh;
    
    
    case lower('activation')
    %=======================================================================
    % [F, winwid, winh] = hewma_gui('activation',index value of image)
    %-----------------------------------------------------------------------
    % Display activation clusters
    
    pop1 = findobj('Tag','ImagePop'); indx = get(pop1,'Value');
    
    % image to define clusters
    imgs = {'hewma_sig.img' 'hewma_sigdiff.img' 'hewma_cp.img' 'gausslongest.img'};
    
    % image to get data
    img2 = {'hewma_t.img' 'hewma_tdiff.img' '' ''};
    
    % check for images
    for i = 1:length(imgs)
        if ~exist(imgs{i}, 'file')
            disp(['Image ' imgs{i} ' cannot be found in current directory.  Please go to a valid HEWMA 2nd-level results directory.']);
            error('Exiting.')
        end
    end
    
    if isempty(img2{indx})
        cl = mask2clusters(imgs{indx}); cluster_orthviews(cl);
    else
        cl = mask2clusters(imgs{indx},img2{indx}); cluster_orthviews(cl,'bivalent');
    
    end
    
    varargout{1} = cl;
    
    % set up timeseries viewer, if asked for
    
    chk = findobj('Tag','TimeseriesChk');
    
    if get(chk,'Value')
        
        % used in button-up fcn callback
        E2 = EXPT;
        clear EXPT

        global VOL
        global f
        global f2
        global EXPT
        EXPT = E2;

        % get coordinate mapping matrix
        VOL = struct('M',cl(1).M);

        % prepare figure
        f1 = figure('Color','w','Name','Hewma plots', 'Tag', 'hewma2display');
        f = f1;     % button callback uses figure f
        
        spmfig = findobj('Tag','Graphics');
        set(spmfig,'WindowButtonUpFcn','[dat, stats, mycov] = hewma_plot_coord_btnupfcn;')
        disp('Timeseries viewer ready.  Click on figure...results returned in dat');
    end


    
    
    
    case lower('classify')
    %=======================================================================
    % [F, winwid, winh] = hewma_gui('classify')
    %-----------------------------------------------------------------------
    % Classify activated voxels
    
    pop1 = findobj('Tag','ClassifyPop'); indx = get(pop1,'Value');
    
    % image to define clusters -- not used yet
    %imgs = {'hewma_sig.img' 'hewma_sigdiff.img' 'gausslongest.img'};
    
    % callback function
    callbk = {'cl = hewma_plot_cpmap;' 'cl = hewma_plot_bivariate;'};
    
    eval(callbk{indx});
    
    varargout{1} = cl;
    
    otherwise
    %=======================================================================
    error('Unknown action string')

    %======================================================================
    
end


return


