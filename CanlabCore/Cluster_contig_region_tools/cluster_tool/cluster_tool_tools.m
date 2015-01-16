function varargout = cluster_tool_tools(varargin)
% CLUSTER_TOOL_TOOLS M-file for cluster_tool_tools.fig
%      CLUSTER_TOOL_TOOLS, by itself, creates a new CLUSTER_TOOL_TOOLS or raises the existing
%      singleton*.
%
%      H = CLUSTER_TOOL_TOOLS returns the handle to a new CLUSTER_TOOL_TOOLS or the handle to
%      the existing singleton*.
%
%      CLUSTER_TOOL_TOOLS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CLUSTER_TOOL_TOOLS.M with the given input arguments.
%
%      CLUSTER_TOOL_TOOLS('Property','Value',...) creates a new CLUSTER_TOOL_TOOLS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cluster_tool_tools_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cluster_tool_tools_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cluster_tool_tools

% Last Modified by GUIDE v2.5 07-Sep-2006 18:11:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cluster_tool_tools_OpeningFcn, ...
                   'gui_OutputFcn',  @cluster_tool_tools_OutputFcn, ...
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


% --- Executes just before cluster_tool_tools is made visible.
function cluster_tool_tools_OpeningFcn(hObject, eventdata, handles, varargin)
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
% varargin   command line arguments to cluster_tool_tools (see VARARGIN)

% Choose default command line output for cluster_tool_tools
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes cluster_tool_tools wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = cluster_tool_tools_OutputFcn(hObject, eventdata, handles) 
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
if isfield(handles,'cl'),cluster_tool('cl',handles.cl,'fdir',handles.fdir);
else cluster_tool('fdir',handles.fdir);end
delete(handles.figure1);


% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
cluster_tool_savefunc


% --- Executes on button press in mask.
function mask_Callback(hObject, eventdata, handles)
cl_mask=cluster_tool_callgui('cluster_tool_masktypeq');
if cl_mask
    handles.pwd=pwd;cd(handles.fdir);
    [FileName,PathName]=uigetfile('*.mat','Select .mat file containing cl structure');
    cd(handles.pwd);
    if ~PathName,return,end
    handles.fdir=PathName;
    load([PathName FileName]);
    try mask=cl;catch,beep,disp('Error: No cl structure in .mat file'),return,end
else
    handles.pwd=pwd;cd(handles.fdir);
    [FileName,PathName]=uigetfile('*.img','Select mask .img');
    cd(handles.pwd);
    if ~PathName,return,end
    handles.fdir=PathName;
    mask=iimg_read_img([PathName FileName],1);
end
if ~isfield(mask,'M')
    if ~isfield(mask,'mat'),beep,disp('Error: no .M or .mat field in source space.'),return,end
    mask.M=mask.mat;
end
if ~isequal(handles.cl(1).voxSize(:),diag(mask(1).M(1:3,1:3)))
    warndlg('Warning: The mask you have selected is not in the same image space as the cl structure. The cl structure will be transferred into the space of the mask. Note that this may lead to undesirable behavior if the mask space has a higher spatial resolution than the space of the cl structure.')
end
if ~cl_mask
    handles.cl=cl_img_intersect(handles.cl,mask);
else
    handles.cl=cl_cl_intersect(mask,handles.cl,'switch');
end
guidata(hObject,handles);


% --- Executes on button press in write.
function write_Callback(hObject, eventdata, handles)
handles.pwd=pwd;cd(handles.fdir);
[FileName,PathName]=uiputfile('*.img','Save .img file');
cd(handles.pwd);
if ~PathName,return,end
handles.fdir=PathName;
use_values=cluster_tool_callgui('cluster_tool_writemasktypeq');
if use_values
    Z=[];
end
XYZ=[];
for k=1:length(handles.cl)
    if size(handles.cl(k).XYZ,1)~=3
        handles.cl(k).XYZ=handles.cl(k).XYZ';
    end
    XYZ(:,end+1:end+handles.cl(k).numVox)=handles.cl(k).XYZ;
    if use_values
        Z(end+1:end+handles.cl(k).numVox)=handles.cl(k).Z;
    end
end
handles.cl(1).nvox=prod(handles.cl(1).dim(1:3));
ind=sub2ind(handles.cl(1).dim(1:3),XYZ(1,:),XYZ(2,:),XYZ(3,:));
dat=zeros(handles.cl(1).nvox,1);
if use_values
    dat(ind)=Z;
else
    dat(ind)=1;
end
if ~isfield(handles.cl,'mat')
    if ~isfield(handles.cl,'M')
        error('No affine matrix in cl structure. Please use the ''Check structure integrity'' tool (cl_check.m).')
    else
        handles.cl(1).mat=handles.cl(1).M;
    end
end
iimg_write_images(dat,handles.cl(1),[PathName FileName]);
% guidata(hObject,handles);

% --- Executes on button press in subdivide.
function subdivide_Callback(hObject, eventdata, handles)
cluster_tool_subdivide('cl',handles.cl,'fdir',handles.fdir);
close(handles.figure1);
% --- Executes on button press in recluster.
function recluster_Callback(hObject, eventdata, handles)


% --- Executes on button press in table.
function table_Callback(hObject, eventdata, handles)
cluster_tool_table('cl',handles.cl,'fdir',handles.fdir);
close(handles.figure1);

% --- Executes on button press in edit.
function edit_Callback(hObject, eventdata, handles)
cluster_tool_edit('cl',handles.cl,'fdir',handles.fdir);
close(handles.figure1);

% --- Executes on button press in check.
function check_Callback(hObject, eventdata, handles)
warndlg('Please watch the matlab command window and answer any questions asked there')
handles.cl=check_cl(handles.cl);
guidata(hObject,handles);



% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
delete(hObject);




% --- Executes on button press in get_betas.
function get_betas_Callback(hObject, eventdata, handles)
cluster_tool_getbetas('cl',handles.cl,'fdir',handles.fdir);
close(handles.figure1);


