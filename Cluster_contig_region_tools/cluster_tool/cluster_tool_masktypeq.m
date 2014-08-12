function varargout = cluster_tool_masktypeq(varargin)
%CLUSTER_TOOL_MASKTYPEQ M-file for cluster_tool_masktypeq.fig
%      CLUSTER_TOOL_MASKTYPEQ, by itself, creates a new CLUSTER_TOOL_MASKTYPEQ or raises the existing
%      singleton*.
%
%      H = CLUSTER_TOOL_MASKTYPEQ returns the handle to a new CLUSTER_TOOL_MASKTYPEQ or the handle to
%      the existing singleton*.
%
%      CLUSTER_TOOL_MASKTYPEQ('Property','Value',...) creates a new CLUSTER_TOOL_MASKTYPEQ using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to cluster_tool_masktypeq_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      CLUSTER_TOOL_MASKTYPEQ('CALLBACK') and CLUSTER_TOOL_MASKTYPEQ('CALLBACK',hObject,...) call the
%      local function named CALLBACK in CLUSTER_TOOL_MASKTYPEQ.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cluster_tool_masktypeq

% Last Modified by GUIDE v2.5 13-Jul-2006 21:50:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cluster_tool_masktypeq_OpeningFcn, ...
                   'gui_OutputFcn',  @cluster_tool_masktypeq_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before cluster_tool_masktypeq is made visible.
function cluster_tool_masktypeq_OpeningFcn(hObject, eventdata, handles, varargin)
movegui(hObject,'center');
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for cluster_tool_masktypeq
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes cluster_tool_masktypeq wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = cluster_tool_masktypeq_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes on button press in overwrite.
function overwrite_Callback(hObject, eventdata, handles)
handles.data=0;
guidata(hObject,handles)
uiresume(handles.figure1)
% hObject    handle to overwrite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in append.
function append_Callback(hObject, eventdata, handles)
handles.data=1;
guidata(hObject,handles)
uiresume(handles.figure1)
% hObject    handle to append (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


