function varargout = cluster_tool_thresh(varargin)
% CLUSTER_TOOL_THRESH M-file for cluster_tool_thresh.fig
%      CLUSTER_TOOL_THRESH, by itself, creates a new CLUSTER_TOOL_THRESH or raises the existing
%      singleton*.
%
%      H = CLUSTER_TOOL_THRESH returns the handle to a new CLUSTER_TOOL_THRESH or the handle to
%      the existing singleton*.
%
%      CLUSTER_TOOL_THRESH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CLUSTER_TOOL_THRESH.M with the given input arguments.
%
%      CLUSTER_TOOL_THRESH('Property','Value',...) creates a new CLUSTER_TOOL_THRESH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cluster_tool_thresh_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cluster_tool_thresh_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cluster_tool_thresh

% Last Modified by GUIDE v2.5 20-Jul-2006 15:43:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cluster_tool_thresh_OpeningFcn, ...
                   'gui_OutputFcn',  @cluster_tool_thresh_OutputFcn, ...
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


% --- Executes just before cluster_tool_thresh is made visible.
function cluster_tool_thresh_OpeningFcn(hObject, eventdata, handles, varargin)
movegui(hObject,'center');
set(handles.title,'string',cluster_tool_getver);
set(handles.t,'value',1);set(handles.p,'value',0);set(handles.dat,'value',0);
set(handles.t_thresh,'value',1);set(handles.p_thresh,'value',0);
handles.thrshtyp='t';handles.imtype='t';
handles.pwd=pwd;handles.fdir=pwd;
for k=1:size(varargin(:),1)
    if strcmp(varargin{k},'cl'),handles.cl=varargin{k+1};
    elseif strcmp(varargin{k},'fdir'),handles.fdir=varargin{k+1};
    elseif strcmp(varargin{k},'img'),handles.img=varargin{k+1};
        set(handles.imagename,'string',['Image: ' handles.img]);
    end
end
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cluster_tool_thresh (see VARARGIN)

% Choose default command line output for cluster_tool_thresh
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes cluster_tool_thresh wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = cluster_tool_thresh_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




% --------------------------------------------------------------------
function imgtype_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to imgtype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.imtype=get(hObject,'string');
if strcmp(handles.imtype,'data')||strcmp(handles.imtype,'p')
    set(handles.threshtype,'visible','off');set(handles.dftxt,'visible','off');set(handles.df,'visible','off');
else set(handles.threshtype,'visible','on');set(handles.dftxt,'visible','on');set(handles.df,'visible','on');
end
if strcmp(handles.imtype,'p'),set(handles.posneg,'visible','off','value',1);
elseif ~strcmp(handles.thrshtyp,'p'),set(handles.posneg,'visible','on');
end
if ~strcmp(handles.imtype,'p'),set(handles.posneg,'visible','on');end
if strcmp(handles.thrshtyp,'t'),handles.thrshtyp='p';set(handles.p_thresh,'value',1);set(handles.t_thresh,'value',0);end
guidata(handles.figure1,handles);



function thresh_Callback(hObject, eventdata, handles)
try handles.threshold=str2double(get(hObject,'string'));catch set(hObject,'string','');handles.threshold=[];end
guidata(hObject,handles);
% hObject    handle to thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresh as text
%        str2double(get(hObject,'String')) returns contents of thresh as a double


% --- Executes during object creation, after setting all properties.
function thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --------------------------------------------------------------------
function threshtype_SelectionChangeFcn(hObject, eventdata, handles)
if strcmp(handles.imtype,'p')
    set(handles.p_thresh,'value',1);set(handles.t_thresh,'value',0);set(handles.posneg,'visible','off','value',1);
else handles.thrshtyp=get(hObject,'string');if ~strcmp(handles.imtype,'t'),set(handles.posneg,'visible','on');end
end
guidata(handles.figure1,handles);
% hObject    handle to threshtype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





function df_Callback(hObject, eventdata, handles)
handles.dfd=str2double(get(hObject,'string'));
if isnan(handles.dfd),set(hObject,'string','');end
guidata(hObject,handles);
% hObject    handle to df (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of df as text
%        str2double(get(hObject,'String')) returns contents of df as a double


% --- Executes during object creation, after setting all properties.
function df_CreateFcn(hObject, eventdata, handles)
% hObject    handle to df (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function k_Callback(hObject, eventdata, handles)
handles.kd=str2double(get(hObject,'string'));
if isnan(handles.kd),set(hObject,'string','');handles.kd=0;end
guidata(hObject,handles);
% hObject    handle to k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k as text
%        str2double(get(hObject,'String')) returns contents of k as a double


% --- Executes during object creation, after setting all properties.
function k_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in imgsel.
function imgsel_Callback(hObject, eventdata, handles)
if isfield(handles,'cl'),handles.append=cluster_tool_callgui('cluster_tool_append');else handles.append=0;end
if isfield(handles,'img'),msgbox('Selected image will be used instead of previously selected image','','warn');end
handles.pwd=pwd;cd(handles.fdir);
[FileName,PathName]=uigetfile('*.img','Select image file.');
cd(handles.pwd);
if isempty(PathName)||isnumeric(PathName),return,end
handles.fdir=PathName;
set(handles.imagename,'string',['Image: ' PathName FileName]);
handles.img=[PathName FileName];
guidata(hObject,handles);
% hObject    handle to imgsel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in done.
function done_Callback(hObject, eventdata, handles)
s='[dat,volInfo,cl]=iimg_threshold(handles.img,''imgtype'',handles.imtype,''thr'',handles.threshold';
if get(handles.posneg,'value')==2&&~strcmp(handles.thrshtyp,'p'),handles.threshold=[-Inf -handles.threshold];
elseif get(handles.posneg,'value')==3,s=[s ',''abs'''];
end
% try
    if ~strcmp(handles.imtype,'data')
        s=[s ',''threshtype'',handles.thrshtyp'];
        if strcmp(handles.thrshtyp,'p')&&strcmp(handles.imtype,'t'),s=[s ',''df'',handles.dfd'];end
    end
    s=[s ',''thr'',handles.threshold'];
    if handles.kd,s=[s ',''k'',handles.kd'];end
    if get(handles.write,'value'),[FileName,PathName]=uiputfile('*.img','Select .img file to write');
        s=[s ',''outnames'',[PathName FileName]'];
    end
    s=[s ');'];
    eval(s)
% catch beep,disp('Error: Not all inputs were specified, please try again.'),return,end
if isfield(handles,'cl'),handles.cl=[handles.cl(:)' cl(:)'];else handles.cl=cl;end
guidata(hObject,handles);
answer=questdlg('Clusters obtained! Return to Acquire/Create Cluster screen?');
if strcmp(answer,'Yes'),ret_Callback(handles.ret,eventdata,handles);end
        
% hObject    handle to done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in ret.
function ret_Callback(hObject, eventdata, handles)
if isfield(handles,'cl'),cluster_tool('cl',handles.cl,'fdir',handles.fdir);
else cluster_tool;end
close(handles.figure1);
    
% hObject    handle to ret (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in write.
function write_Callback(hObject, eventdata, handles)
% hObject    handle to write (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of write




% --- Executes on selection change in posneg.
function posneg_Callback(hObject, eventdata, handles)
% hObject    handle to posneg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns posneg contents as cell array
%        contents{get(hObject,'Value')} returns selected item from posneg


% --- Executes during object creation, after setting all properties.
function posneg_CreateFcn(hObject, eventdata, handles)
set(hObject,'string','Positive|Negative|Both');
% hObject    handle to posneg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function imagename_CreateFcn(hObject,eventdata,handles)
set(hObject,'string',['Image: ' handles.img]);
