function varargout = cluster_tool_getbetas(varargin)
% CLUSTER_TOOL_GETBETAS M-file for cluster_tool_getbetas.fig
%      CLUSTER_TOOL_GETBETAS, by itself, creates a new CLUSTER_TOOL_GETBETAS or raises the existing
%      singleton*.
%
%      H = CLUSTER_TOOL_GETBETAS returns the handle to a new CLUSTER_TOOL_GETBETAS or the handle to
%      the existing singleton*.
%
%      CLUSTER_TOOL_GETBETAS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CLUSTER_TOOL_GETBETAS.M with the given input arguments.
%
%      CLUSTER_TOOL_GETBETAS('Property','Value',...) creates a new CLUSTER_TOOL_GETBETAS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cluster_tool_getbetas_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cluster_tool_getbetas_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cluster_tool_getbetas

% Last Modified by GUIDE v2.5 07-Sep-2006 19:22:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cluster_tool_getbetas_OpeningFcn, ...
                   'gui_OutputFcn',  @cluster_tool_getbetas_OutputFcn, ...
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


% --- Executes just before cluster_tool_getbetas is made visible.
function cluster_tool_getbetas_OpeningFcn(hObject, eventdata, handles, varargin)
movegui(hObject,'center');
set(handles.title,'string',cluster_tool_getver);
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
% varargin   command line arguments to cluster_tool_getbetas (see VARARGIN)

% Choose default command line output for cluster_tool_getbetas
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes cluster_tool_getbetas wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = cluster_tool_getbetas_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in expt.
function expt_Callback(hObject, eventdata, handles)
try
    handles.cl=extract_contrast_data(evalin('base','EXPT.SNPM.P'),handles.cl);
catch beep,disp('Error: Could not get EXPT.SNPM.P from base workspace')
end
guidata(hObject,handles);


% --- Executes on button press in get_from_img.
function get_from_img_Callback(hObject, eventdata, handles)
try p = spm_get([0 Inf],'*.img','Select .img(s) to extract data from',handles.fdir,0);catch end
if ~isempty(p)
    for k=1:size(p,1)
        P{k}=deblank(p(k,:));
    end
    handles.cl=extract_contrast_data(P,handles.cl);
end
guidata(hObject,handles);

% --- Executes on button press in quit.
function quit_Callback(hObject, eventdata, handles)
delete(handles.figure1);


% --- Executes on button press in ret.
function ret_Callback(hObject, eventdata, handles)
if isfield(handles,'cl'),cluster_tool_tools('cl',handles.cl,'fdir',handles.fdir);
else cluster_tool_tools('fdir',handles.fdir);end
delete(handles.figure1);


% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
cluster_tool_savefunc


