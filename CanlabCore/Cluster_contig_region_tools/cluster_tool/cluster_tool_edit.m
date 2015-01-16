function varargout = cluster_tool_edit(varargin)
% CLUSTER_TOOL_EDIT M-file for cluster_tool_edit.fig
%      CLUSTER_TOOL_EDIT, by itself, creates a new CLUSTER_TOOL_EDIT or raises the existing
%      singleton*.
%
%      H = CLUSTER_TOOL_EDIT returns the handle to a new CLUSTER_TOOL_EDIT or the handle to
%      the existing singleton*.
%
%      CLUSTER_TOOL_EDIT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CLUSTER_TOOL_EDIT.M with the given input arguments.
%
%      CLUSTER_TOOL_EDIT('Property','Value',...) creates a new CLUSTER_TOOL_EDIT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cluster_tool_edit_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cluster_tool_edit_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cluster_tool_edit

% Last Modified by GUIDE v2.5 07-Sep-2006 17:34:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cluster_tool_edit_OpeningFcn, ...
                   'gui_OutputFcn',  @cluster_tool_edit_OutputFcn, ...
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


% --- Executes just before cluster_tool_edit is made visible.
function cluster_tool_edit_OpeningFcn(hObject, eventdata, handles, varargin)
movegui(hObject,'center');
set(handles.title,'string',[cluster_tool_getver]);
handles.pwd=pwd;handles.fdir=pwd;
for k=1:size(varargin(:),1)
    if strcmp(varargin{k},'cl'),handles.cl=varargin{k+1};
    elseif strcmp(varargin{k},'fdir'),handles.fdir=varargin{k+1};
    end
end

s='';
for k=1:length(handles.cl)
    if k~=1
        s=[s '|'];
    end
    s=[s 'Cluster of ' num2str(handles.cl(k).numVox) ' voxels with center ' num2str(handles.cl(k).mm_center(1)) ', ' num2str(handles.cl(k).mm_center(2)) ', ' num2str(handles.cl(k).mm_center(3)) '.'];
end
set(handles.cl_list,'String',s);set(handles.cl_list,'Max',length(handles.cl));
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cluster_tool_edit (see VARARGIN)

% Choose default command line output for cluster_tool_edit
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes cluster_tool_edit wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = cluster_tool_edit_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in cl_list.
function cl_list_Callback(hObject, eventdata, handles)
% hObject    handle to cl_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns cl_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cl_list


% --- Executes during object creation, after setting all properties.
function cl_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cl_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




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


% --- Executes on button press in del.
function del_Callback(hObject, eventdata, handles)
handles.cl(get(handles.cl_list,'Value'))=[];
set(handles.cl_list,'Max',length(handles.cl));
s=get(handles.cl_list,'String');
s(get(handles.cl_list,'Value'),:)=[];
set(handles.cl_list,'Value',1);
set(handles.cl_list,'String',s);
guidata(hObject,handles);

