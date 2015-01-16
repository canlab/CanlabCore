function varargout = cluster_tool_subdivide(varargin)
% CLUSTER_TOOL_SUBDIVIDE M-file for cluster_tool_subdivide.fig
%      CLUSTER_TOOL_SUBDIVIDE, by itself, creates a new CLUSTER_TOOL_SUBDIVIDE or raises the existing
%      singleton*.
%
%      H = CLUSTER_TOOL_SUBDIVIDE returns the handle to a new CLUSTER_TOOL_SUBDIVIDE or the handle to
%      the existing singleton*.
%
%      CLUSTER_TOOL_SUBDIVIDE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CLUSTER_TOOL_SUBDIVIDE.M with the given input arguments.
%
%      CLUSTER_TOOL_SUBDIVIDE('Property','Value',...) creates a new CLUSTER_TOOL_SUBDIVIDE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cluster_tool_subdivide_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cluster_tool_subdivide_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cluster_tool_subdivide

% Last Modified by GUIDE v2.5 04-Sep-2006 11:37:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cluster_tool_subdivide_OpeningFcn, ...
                   'gui_OutputFcn',  @cluster_tool_subdivide_OutputFcn, ...
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


% --- Executes just before cluster_tool_subdivide is made visible.
function cluster_tool_subdivide_OpeningFcn(hObject, eventdata, handles, varargin)
movegui(hObject,'center');
set(handles.title,'string',[cluster_tool_getver]);
handles.pwd=pwd;handles.fdir=pwd;
for k=1:size(varargin(:),1)
    if strcmp(varargin{k},'cl'),handles.cl=varargin{k+1};
    elseif strcmp(varargin{k},'fdir'),handles.fdir=varargin{k+1};
    end
end
handles.criterion=0;handles.k=1;handles.min_adj=0;
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cluster_tool_subdivide (see VARARGIN)

% Choose default command line output for cluster_tool_subdivide
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes cluster_tool_subdivide wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = cluster_tool_subdivide_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


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


% --- Executes on selection change in adjacency.
function adjacency_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function adjacency_CreateFcn(hObject, eventdata, handles)
set(hObject,'String','Surface|Edge|Corner');
set(hObject,'Value',2);
disp('**If a warning prints, you can safely ignore it.**')
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in done.
function done_Callback(hObject, eventdata, handles)
if get(handles.adjacency,'Value')==1
    adjacency='surface';
elseif get(handles.adjacency,'Value')==2
    adjacency='edge';
elseif get(handles.adjacency,'Value')==3
    adjacency='corner';
end
cl=cl_subdivide(handles.cl,handles.criterion,handles.min_adj,adjacency,handles.k);
if ~isempty(cl)
    handles.cl=cl;
else
    warndlg('Resulting ''cl'' structure was empty, the initial structure is being kept in memory.');
end
guidata(hObject,handles);

function criterion_Callback(hObject, eventdata, handles)
try handles.criterion=str2double(get(hObject,'string'));catch set(hObject,'string','0');handles.criterion=0;end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function criterion_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function k_Callback(hObject, eventdata, handles)
try handles.k=str2double(get(hObject,'string'));catch set(hObject,'string','1');handles.k=1;end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function k_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function minadj_Callback(hObject, eventdata, handles)
try handles.min_adj=str2double(get(hObject,'string'));catch set(hObject,'string','0');handles.min_adj=0;end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function minadj_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minadj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


