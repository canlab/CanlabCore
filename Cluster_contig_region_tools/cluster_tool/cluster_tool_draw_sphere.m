function varargout = cluster_tool_draw_sphere(varargin)
% CLUSTER_SRC_TOOL_DRAW_SPHERE M-file for cluster_src_tool_draw_sphere.fig
%      CLUSTER_SRC_TOOL_DRAW_SPHERE, by itself, creates a new CLUSTER_SRC_TOOL_DRAW_SPHERE or raises the existing
%      singleton*.
%
%      H = CLUSTER_SRC_TOOL_DRAW_SPHERE returns the handle to a new CLUSTER_SRC_TOOL_DRAW_SPHERE or the handle to
%      the existing singleton*.
%
%      CLUSTER_SRC_TOOL_DRAW_SPHERE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CLUSTER_SRC_TOOL_DRAW_SPHERE.M with the given input arguments.
%
%      CLUSTER_SRC_TOOL_DRAW_SPHERE('Property','Value',...) creates a new CLUSTER_SRC_TOOL_DRAW_SPHERE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cluster_src_tool_draw_sphere_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cluster_src_tool_draw_sphere_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cluster_src_tool_draw_sphere

% Last Modified by GUIDE v2.5 01-Aug-2006 16:11:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cluster_tool_draw_sphere_OpeningFcn, ...
                   'gui_OutputFcn',  @cluster_tool_draw_sphere_OutputFcn, ...
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


% --- Executes just before cluster_src_tool_draw_sphere is made visible.
function cluster_tool_draw_sphere_OpeningFcn(hObject, eventdata, handles, varargin)
movegui(hObject,'center');
set(handles.title,'string',cluster_tool_getver);
handles.src='Image';
handles.intype='Manual';
handles.pwd=pwd;handles.fdir=pwd;
for k=1:size(varargin(:),1)
    if strcmp(varargin{k},'cl'),handles.cl=varargin{k+1};
    elseif strcmp(varargin{k},'fdir'),handles.fdir=varargin{k+1};
    elseif strcmp(varargin{k},'img'),handles.img=varargin{k+1};
        set(handles.sourcespace,'string',['Image: ' handles.img]);
    end
end
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cluster_src_tool_draw_sphere (see VARARGIN)

% Choose default command line output for cluster_src_tool_draw_sphere
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes cluster_src_tool_draw_sphere wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = cluster_tool_draw_sphere_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function radius_Callback(hObject, eventdata, handles)
try handles.r=str2double(get(hObject,'string'));catch set(hObject,'string','');handles.r=[];end
guidata(hObject,handles);
% Hints: get(hObject,'String') returns contents of radius as text
%        str2double(get(hObject,'String')) returns contents of radius as a double


% --- Executes during object creation, after setting all properties.
function radius_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in selector.
function selector_Callback(hObject, eventdata, handles)
if strcmp(handles.src,'Image')
    if isfield(handles,'img'),msgbox('Selected image will be used instead of previously selected image','','warn');end
    handles.pwd=pwd;cd(handles.fdir);
    [FileName,PathName]=uigetfile('*.img','Select image file.');
    cd(handles.pwd);
    if ~PathName,return,end
    handles.fdir=PathName;
    set(handles.sourcespace,'String',['Source: ' PathName FileName]);
    handles.img=[PathName FileName];
else
    handles.pwd=pwd;cd(handles.fdir);
    [FileName,PathName]=uigetfile('*.mat','Select mat file.');
    cd(handles.pwd);
    if ~PathName,return,end
    handles.fdir=PathName;
    load([PathName FileName]);
    try 
        handles.cl_space=cl;
        set(handles.sourcespace,'String',['Source: ' PathName FileName]);
        handles.cl_space(1).fname=[PathName FileName];
    catch
        beep,disp('Error: .mat file likely does not contain a cl structure'),return
    end
end
guidata(hObject,handles);



% --- Executes on button press in done.
function done_Callback(hObject, eventdata, handles)
try
    if get(handles.image_src,'Value')
        cl_space=iimg_read_img(handles.img,1);
    else
        cl_space=handles.cl_space(1);
    end
    if ~isfield(cl_space,'M')
        if ~isfield(cl_space,'mat'),beep,disp('Error: no .M or .mat field in source space.'),return,end
        M=cl_space.mat;
    else
        M=cl_space.M;
    end
    if get(handles.Man_in,'Value')
        XYZmm=sphere_3d([handles.x_val handles.y_val handles.z_val]',handles.r,diag(M(1:3,1:3))');
    else
        XYZmm=sphere_3d(handles.sphere_centers,handles.r,diag(M(1:3,1:3))');
    end
    XYZ=mmToVoxel(XYZmm,M,'valid');
    cl_ind=spm_clusters(XYZ);
    for k=1:max(cl_ind)
        cl(k).XYZ=XYZ(:,find(cl_ind==k));
        cl(k).XYZmm=voxelToMm(cl(k).XYZ,M);
        cl(k).numVox=length(find(cl_ind==k));
        cl(k).M=M;
        cl(k).voxSize=diag(cl(k).M(1:3,1:3));
        cl(k).name='';
        cl(k).Z(1:cl(k).numVox)=1;
        cl(k).mm_center=center_of_mass(cl(k).XYZmm,cl(k).Z);
        cl(k).dim=cl_space.dim;
        cl(k).threshold=[1 NaN];
        cl(k).title=['Cluster of ' num2str(cl(k).numVox) ' voxels, created by cluster_tool_draw_sphere GUI. Image space from: ' cl_space.fname];
    end
    if get(handles.domask,'Value')
        if get(handles.image_src,'Value')
            cl=cl_img_intersect(cl,handles.img);
        else
            cl=cl_cl_intersect(cl,handles.cl_space);
        end
    end     
    if isfield(handles,'cl')
        if isfield(handles,'cl'),append=cluster_tool_callgui('cluster_tool_append');else handles.append=0;end
        if append
            handles.cl=[handles.cl(:);cl(:)];
        else
            handles.cl=cl;
        end
    else
        handles.cl=cl;
    end
catch beep,disp('Error: Not all required data is specified. Check that you have selected a file, specified a radius, and specified valid sphere centers.'),return
end
guidata(handles.figure1,handles);


% --- Executes on button press in ret.
function ret_Callback(hObject, eventdata, handles)
if isfield(handles,'cl'),cluster_tool('cl',handles.cl,'fdir',handles.fdir);
else cluster_tool('fdir',handles.fdir);end
close(handles.figure1);
    

function x_Callback(hObject, eventdata, handles)
try handles.x_val=str2double(get(hObject,'string'));catch set(hObject,'string','');handles.x_val=[];end
guidata(hObject,handles);
% Hints: get(hObject,'String') returns contents of x as text
%        str2double(get(hObject,'String')) returns contents of x as a double


% --- Executes during object creation, after setting all properties.
function x_CreateFcn(hObject, eventdata, handles)


% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function y_Callback(hObject, eventdata, handles)
try handles.y_val=str2double(get(hObject,'string'));catch set(hObject,'string','');handles.y_val=[];end
guidata(hObject,handles);

% Hints: get(hObject,'String') returns contents of y as text
%        str2double(get(hObject,'String')) returns contents of y as a double


% --- Executes during object creation, after setting all properties.
function y_CreateFcn(hObject, eventdata, handles)


% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function z_Callback(hObject, eventdata, handles)
try handles.z_val=str2double(get(hObject,'string'));catch set(hObject,'string','');handles.z_val=[];end
guidata(hObject,handles);

% Hints: get(hObject,'String') returns contents of z as text
%        str2double(get(hObject,'String')) returns contents of z as a double


% --- Executes during object creation, after setting all properties.
function z_CreateFcn(hObject, eventdata, handles)


% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function var_name_Callback(hObject, eventdata, handles)
try 
    handles.sphere_centers=evalin('base',get(hObject,'String'));
catch beep,disp('Error: Could not get data from workspace variable'),return,end
if size(size(handles.sphere_centers))~=2||~isnumeric(handles.sphere_centers)||(min(size(handles.sphere_centers))~=3&max(size(handles.sphere_centers))~=3)
    beep,disp('Error: Specified variable does not appear to be a valid n x 3 or 3 x n array of center coordinates'),return
end
if size(handles.sphere_centers,1)~=3
    handles.sphere_centers=handles.sphere_centers';
elseif size(handles.sphere_centers,2)==3
    beep,disp('Warning: Assuming that X Y Z values are in ROWS, not COLUMNS (that is, that each COLUMN is a point).'),disp('If this is not the case, transpose the matrix and try again')
end
guidata(hObject,handles);

% Hints: get(hObject,'String') returns contents of var_name as text
%        str2double(get(hObject,'String')) returns contents of var_name as a double


% --- Executes during object creation, after setting all properties.
function var_name_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --------------------------------------------------------------------
function source_selector_SelectionChangeFcn(hObject, eventdata, handles)
handles.src=get(hObject,'String');
if strcmp(handles.src,'Image')
    set(handles.selector,'String','Select Image');
    if isfield(handles,'img')
        set(handles.sourcespace,'String',['Source: ' handles.img]);
    else
        set(handles.sourcespace,'string','Source: (none selected)');
    end
else
    set(handles.selector,'String','Select File');
    if isfield(handles,'cl_space')
        set(handles.sourcespace,'String',['Source: ' handles.cl_space.fname]);
    else
        set(handles.sourcespace,'String','Source: (none selected)');
    end
end
guidata(handles.figure1,handles);

% --------------------------------------------------------------------
function center_input_selector_SelectionChangeFcn(hObject, eventdata, handles)
handles.intype=get(hObject,'String');
if strcmp(handles.intype,'Manual')
    set(handles.x,'Visible','on');
    set(handles.y,'Visible','on');
    set(handles.z,'Visible','on');
    set(handles.x_text,'Visible','on');
    set(handles.y_text,'Visible','on');
    set(handles.z_text,'Visible','on');
    set(handles.var_name,'Visible','off');
    set(handles.center_text,'String','Sphere Center:');
    if isfield(handles,'wrkspc_var'),handles=rmfield(handles,'wrkspc_var');end
else
    set(handles.x,'Visible','off');
    set(handles.y,'Visible','off');
    set(handles.z,'Visible','off');
    set(handles.x_text,'Visible','off');
    set(handles.y_text,'Visible','off');
    set(handles.z_text,'Visible','off');
    set(handles.var_name,'Visible','on');
    set(handles.center_text,'String','Variable Name:');
    if isfield(handles,'x_val'),handles=rmfield(handles,'x_val');end
    if isfield(handles,'y_val'),handles=rmfield(handles,'y_val');end
    if isfield(handles,'z_val'),handles=rmfield(handles,'z_val');end
end
guidata(handles.figure1,handles);




