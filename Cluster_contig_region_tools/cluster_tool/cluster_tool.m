function varargout = cluster_tool(varargin)
% CLUSTER_TOOL M-file for cluster_tool.fig
%      CLUSTER_TOOL, by itself, creates a new CLUSTER_TOOL or raises the existing
%      singleton*.
%
%      H = CLUSTER_TOOL returns the handle to a new CLUSTER_TOOL or the handle to
%      the existing singleton*.
%
%      CLUSTER_TOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CLUSTER_TOOL.M with the given input arguments.
%
%      CLUSTER_TOOL('Property','Value',...) creates a new CLUSTER_TOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cluster_tool_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cluster_tool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cluster_tool

% Last Modified by GUIDE v2.5 07-Sep-2006 19:41:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cluster_tool_OpeningFcn, ...
                   'gui_OutputFcn',  @cluster_tool_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before cluster_tool is made visible.
function cluster_tool_OpeningFcn(hObject, eventdata, handles, varargin)
movegui(hObject,'center');
set(handles.title,'string',[cluster_tool_getver]);
handles.pwd=pwd;handles.fdir=pwd;
for k=1:size(varargin(:),1)
    if strcmp(varargin{k},'cl'),handles.cl=varargin{k+1};
    elseif strcmp(varargin{k},'fdir'),handles.fdir=varargin{k+1};
    end
end


% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cluster_tool (see VARARGIN)

% Choose default command line output for cluster_tool
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes cluster_tool wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = cluster_tool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in open_image.
function open_image_Callback(hObject, eventdata, handles)
if isfield(handles,'cl'),handles.append=cluster_tool_callgui('cluster_tool_append');else handles.append=0;end
handles.pwd=pwd;cd(handles.fdir);
[FileName,PathName]=uigetfile('*.img','Select image file.');
cd(handles.pwd);
if isempty(FileName)||~ischar(FileName),beep,disp('No file found'),return,end
handles.fdir=PathName;
if handles.append,cluster_tool_thresh('img',[PathName FileName],'fdir',handles.fdir,'cl',handles.cl);
else cluster_tool_thresh('img',[PathName FileName],'fdir',handles.fdir);
end
guidata(handles.figure1,handles);
close(handles.figure1);


% --- Executes on button press in get_mat.
function get_mat_Callback(hObject, eventdata, handles)
if isfield(handles,'cl'),handles.append=cluster_tool_callgui('cluster_tool_append');else handles.append=0;end
handles.pwd=pwd;cd(handles.fdir);
[FileName,PathName]=uigetfile('*.mat','Select .mat file containing cl structure');
cd(handles.pwd);
if ~PathName,return,end
handles.fdir=PathName;
load([PathName FileName]);
if ~exist('cl','var')||~isstruct(cl),beep,disp('Error: .mat file does not contain ''cl'' structure.'),guidata(hObject,handles);return,end
if handles.append,handles.cl=[handles.cl(:)' cl(:)'];else handles.cl=cl;end
guidata(hObject,handles);


% --- Executes on button press in quit.
function quit_Callback(hObject, eventdata, handles)

close;


% --- Executes on button press in tools.
function tools_Callback(hObject, eventdata, handles)

try 
    cluster_tool_tools('cl',handles.cl,'fdir',handles.fdir);
catch
    beep,disp('Error: No ''cl'' structure in memory.'),return,
end
close(handles.figure1);


% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
cluster_tool_savefunc


% --- Executes on button press in mask.
function mask_Callback(hObject, eventdata, handles)
if isfield(handles,'cl'),handles.append=cluster_tool_callgui('cluster_tool_append');else handles.append=0;end
mask=cluster_tool_callgui('cluster_tool_maskq');
if mask
    handles.pwd=pwd;cd(handles.fdir);
    [FileName,PathName]=uigetfile('*.img','Select .img with DATA to extract...');
    cd(handles.pwd);if isempty(PathName),return,end
    handles.fdir=PathName;
    answer=inputdlg('If t image, enter df (otherwise cancel):');
    if isempty(answer),cl=mask2clusters([],[PathName FileName]);
    else cl=mask2clusters([],[PathName FileName],str2double(answer));end
else answer=inputdlg('If t image, enter df (otherwise cancel):');
    if isempty(answer),cl=mask2clusters([]);
    else cl=mask2clusters([],str2double(answer));end
end
if handles.append,handles.cl=[handles.cl(:)' cl(:)'],else handles.cl=cl;end
guidata(hObject,handles);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% Hint: delete(hObject) closes the figure
delete(hObject);




% --- Executes on button press in draw.
function draw_Callback(hObject, eventdata, handles)
if isfield(handles,'cl'),handles.append=cluster_tool_callgui('cluster_tool_append');else handles.append=0;end
if handles.append,cluster_tool_draw_sphere('cl',handles.cl,'fdir',handles.fdir);else cluster_tool_draw_sphere('fdir',handles.fdir);end
close(handles.figure1);




% --- Executes on button press in get_workspc.
function get_workspc_Callback(hObject, eventdata, handles)
varname=inputdlg('Enter variable name');
if ~isempty(varname)
    handles.cl=evalin('base',varname{1});
end
guidata(hObject,handles);


