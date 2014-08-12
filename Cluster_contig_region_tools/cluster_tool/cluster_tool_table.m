function varargout = cluster_tool_table(varargin)
% CLUSTER_TOOL_TABLE M-file for cluster_tool_table.fig
%      CLUSTER_TOOL_TABLE, by itself, creates a new CLUSTER_TOOL_TABLE or raises the existing
%      singleton*.
%
%      H = CLUSTER_TOOL_TABLE returns the handle to a new CLUSTER_TOOL_TABLE or the handle to
%      the existing singleton*.
%
%      CLUSTER_TOOL_TABLE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CLUSTER_TOOL_TABLE.M with the given input arguments.
%
%      CLUSTER_TOOL_TABLE('Property','Value',...) creates a new CLUSTER_TOOL_TABLE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cluster_tool_table_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cluster_tool_table_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cluster_tool_table

% Last Modified by GUIDE v2.5 17-Oct-2006 10:33:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cluster_tool_table_OpeningFcn, ...
                   'gui_OutputFcn',  @cluster_tool_table_OutputFcn, ...
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


% --- Executes just before cluster_tool_table is made visible.
function cluster_tool_table_OpeningFcn(hObject, eventdata, handles, varargin)
movegui(hObject,'center');
set(handles.title,'string',[cluster_tool_getver]);
handles.pwd=pwd;handles.fdir=pwd;
for k=1:size(varargin(:),1)
    if strcmp(varargin{k},'cl'),handles.cl=varargin{k+1};
    elseif strcmp(varargin{k},'fdir'),handles.fdir=varargin{k+1};
    end
end
% s={'Print Labels (Talairach atlas--coordinates converted','from MNI using mni2tal)'};
% set(handles.tal_labels,'String',s);


% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cluster_tool_table (see VARARGIN)

% Choose default command line output for cluster_tool_table
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes cluster_tool_table wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = cluster_tool_table_OutputFcn(hObject, eventdata, handles) 
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


% --- Executes on button press in write.
function write_Callback(hObject, eventdata, handles)


% --- Executes on button press in do_subcl.
function do_subcl_Callback(hObject, eventdata, handles)


% --- Executes on button press in do_labels.
function do_labels_Callback(hObject, eventdata, handles)
if get(handles.do_labels,'Value')
    set(handles.tal_labels,'Value',0)
end

% --- Executes on button press in done.
function done_Callback(hObject, eventdata, handles)
if get(handles.tal_labels,'Value')
    for i=length(handles.cl)
        handles.cl(i).XYZmm=mni2tal(handles.cl(i).XYZmm);
    end
end
s='cluster_table(handles.cl,get(handles.do_subcl,''Value''),any([get(handles.do_labels,''Value'') get(handles.tal_labels,''Value'')])';
if get(handles.write,'Value')
    handles.pwd=pwd;cd(handles.fdir);
    [FileName,PathName]=uiputfile('*.tsv','Save as tab seperated values...');
    cd(handles.pwd);
    if isempty(FileName)||~ischar(FileName),beep,disp('No file found'),return,end
    handles.fdir=PathName;
    s=[s ',''writefile'',[PathName FileName]'];
end
if get(handles.do_labels,'Value')||get(handles.tal_labels,'Value')
    try
        xyz=evalin('base','xyz');
        L3=evalin('base','L3');
        L5=evalin('base','L5');
        s=[s ',''tal_info'',xyz,L3,L5'];
    catch
        if get(handles.tal_labels,'Value')
            s=[s ',''talairach'''];
        end
    end
end
s=[s ');'];
eval(s)
beep,disp('Done printing table.')
  

function tal_labels_CreateFcn(hObject, eventdata,handles)
% s={'Print Labels (Talairach atlas--coordinates converted','from MNI using mni2tal)'};
% set(handles.tal_labels,'String',s)


% --- Executes on button press in tal_labels.
function tal_labels_Callback(hObject, eventdata, handles)
if get(handles.tal_labels,'Value')
    set(handles.do_labels,'Value',0)
end


