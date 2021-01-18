function jframe = getjframe(hFig)
% getjframe Retrieves the underlying Java frame for a figure
%
% Syntax:
%    jframe = getjframe(hFigureOrComponent)
%
% Description:
%    GETJFRAME retrieves the current figure (gcf)'s underlying Java frame,
%    thus enabling access to all 35 figure callbacks that are not exposed
%    by Matlab's figure.
%
%    Notable callbacks include: FocusGainedCallback, FocusLostCallback, 
%    KeyPressedCallback, KeyReleasedCallback, MouseEnteredCallback, 
%    MouseExitedCallback, MousePressedCallback, MouseReleasedCallback,
%    WindowActivatedCallback, WindowClosedCallback, WindowClosingCallback,
%    WindowOpenedCallback, WindowStateChangedCallback and 22 others.
%
%    The returned jframe object also allows access to other useful window
%    features: 'AlwaysOnTop', 'CloseOnEscapeEnabled', 'Resizable',
%    'Enabled', 'HWnd' (for those interested in Windows integration) etc.
%    Type "get(jframe)" to see the full list of properties.
%
%    GETJFRAME(hFig) retrieves a specific figure's underlying Java frame.
%    hFig is a Matlab handle, or a list of graphic handles (not necessarily
%    figure handle(s) - the handles' containing figure is used).
%
% Examples:
%    get(getjframe,'ListOfCallbacks');  %display list of supported callbacks
%    set(getjframe,'WindowStateChangedCallback','disp(''Window min/maxed'')')
%    set(getjframe,'WindowDeiconifiedCallback',@winMaximizedCallback)
%    set(getjframe,'WindowIconifiedCallback',{@winMinimizedCallback,mydata})
%    set(getjframe,'CloseOnEscapeEnabled','on')
%    jframes = getjframe([gcf,hButton]);  % get 2 java frames, from 2 figures
%
% Bugs and suggestions:
%    Please send to Yair Altman (altmany at gmail dot com)
%
% Change log:
%    2007-08-05: First version posted on <a href="http://www.mathworks.com/matlabcentral/fileexchange/loadAuthor.do?objectType=author&mfx=1&objectId=1096533#">MathWorks File Exchange</a>
%    2007-08-11: Added Matlab figure handle property; improved responsiveness handling; added support for array of handles; added sanity checks for illegal handles
%    2011-10-14: Fix for R2011b
%    2014-09-02: Fix for R2014b (HG2)
%    2018-06-28: Additional fixes for R2014b onward (HG2)
%    2019-03-01: Extra check for additional robustness
%
% See also:
%    gcf, findjobj (on the File Exchange)

% Programmed by Yair M. Altman: altmany(at)gmail.com
% $Revision: 1.2 $  $Date: 2020/02/07 07:00:37 $

  try
      % Default figure = current (gcf)
      if nargin < 1 || ~all(ishghandle(hFig))
          if nargin && ~all(ishghandle(hFig))
              error('hFig must be a valid GUI handle or array of handles');
          end
          hFig = gcf;
      end

      % Require Java engine to run
      if ~usejava('jvm')
          error([mfilename ' requires Java to run.']);
      end

      % Initialize output var (needed in case hFig is empty)
      jframe = handle([]);

      % Loop over all requested figures
      for figIdx = 1 : length(hFig)

          % Get the root Java frame
          jff = getJFrame(hFig(figIdx));

          % Add the 'ListOfCallbacks' read-only property to jframe
          fields = fieldnames(get(jff));
          cb1Idx = find(~cellfun('isempty',strfind(fields,'Callback'))); %#ok<*STRCL1>
          cb2Idx = find(~cellfun('isempty',strfind(fields,'CallbackData')));
          cbIdx = setdiff(cb1Idx, cb2Idx)';
          cbNames = fields(cbIdx);
          jframe(figIdx) = handle(jff,'callbackproperties');
          addReadOnlyProp(jframe(figIdx),'ListOfCallbacks',cbNames);

          % Add another read-only property for direct access to the wrapped java object
          addReadOnlyProp(jframe(figIdx),'JavaComponent',jff);

          % Add another read-only property for direct access to the original Matlab figure handle
          addReadOnlyProp(jframe(figIdx),'MatlabFigureHandle',hFig(figIdx));
      end

  % Error handling
  catch
      v = version;
      if v(1)<='6'
          err.message = lasterr;  %#ok<*LERR> % no lasterror function...
      else
          err = lasterror;
      end
      try
          err.message = regexprep(err.message,'Error using ==> [^\n]+\n','');
      catch
          try
              % Another approach, used in Matlab 6 (where regexprep is unavailable)
              startIdx = findstr(err.message,'Error using ==> '); %#ok<*FSTR>
              stopIdx = findstr(err.message,char(10)); %#ok<*CHARTEN>
              for idx = length(startIdx) : -1 : 1
                  idx2 = min(find(stopIdx > startIdx(idx)));  %#ok ML6
                  err.message(startIdx(idx):stopIdx(idx2)) = [];
              end
          catch
              % never mind...
          end
      end
      if isempty(findstr(mfilename,err.message))
          % Indicate error origin, if not already stated within the error message
          err.message = [mfilename ': ' err.message];
      end
      if v(1)<='6'
          while err.message(end)==char(10)
              err.message(end) = [];  % strip excessive Matlab 6 newlines
          end
          error(err.message);
      else
          rethrow(err);
      end
  end

%% Add a read-only property to an object
function addReadOnlyProp(obj,propName,initValue)
  try
      % This will fail in R2014b+ since schema.prop is no longer supported in HG2
      sp = schema.prop(obj,propName,'mxArray');
      set(obj,propName,initValue);
      set(sp,'AccessFlags.PublicSet','off');
  catch
      % Never mind - property might already exist...
  end

%% Get the root Java frame (up to 10 tries, to wait for figure to become responsive)
function jframe = getJFrame(hFigHandle)

  % Ensure that hFig is a figure handle...
  hFig = ancestor(hFigHandle,'figure');
  if isempty(hFig)
      error(['Cannot retrieve the figure handle for handle ' num2str(hFigHandle)]);
  end

  jframe = [];
  maxTries = 10;
  while maxTries > 0
      try
          % Get the figure's underlying Java frame
          oldWarn = warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
          jf = get(handle(hFig),'JavaFrame');
          warning(oldWarn);

          % Get the Java frame's root frame handle
          %jframe = jf.getFigurePanelContainer.getComponent(0).getRootPane.getParent;
          try
              jClient = jf.fFigureClient;  % This works up to R2011a
          catch
              try
                  jClient = jf.fHG1Client;  % This works from R2008b-R2014a
              catch
                  jClient = jf.fHG2Client;  % This works from R2014b and up
              end
          end
          jframe = jClient.getWindow;  % equivalent to above...
          if ~isempty(jframe)
              break;
          else
              maxTries = maxTries - 1;
              drawnow; pause(0.1);
          end
      catch
          maxTries = maxTries - 1;
          drawnow; pause(0.1);
      end
  end
  if isempty(jframe)
      try jframe = javax.swing.SwingUtilities.getWindowAncestor(jf.getAxisComponent); catch, end % see controllibutils.centerfig (%matlabroot%\toolbox\shared\controllib\graphics\+controllibutils.m)
  end
  if isempty(jframe)
      error('Cannot retrieve the java frame for the specified handle');
  end
