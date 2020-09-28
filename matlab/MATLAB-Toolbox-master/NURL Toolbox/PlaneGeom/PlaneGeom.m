function varargout = PlaneGeom(varargin)
% PLANEGEOM MATLAB code for PlaneGeom.fig
%      PLANEGEOM, by itself, creates a new PLANEGEOM or raises the existing
%      singleton*.
% 
%      H = PLANEGEOM returns the handle to a new PLANEGEOM or the handle to
%      the existing singleton*.
%
%      PLANEGEOM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLANEGEOM.M with the given input arguments.
%
%      PLANEGEOM('Property','Value',...) creates a new PLANEGEOM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PlaneGeom_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PlaneGeom_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PlaneGeom

% Last Modified by GUIDE v2.5 03-Sep-2016 15:09:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PlaneGeom_OpeningFcn, ...
                   'gui_OutputFcn',  @PlaneGeom_OutputFcn, ...
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


% --- Executes just before PlaneGeom is made visible.
function PlaneGeom_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PlaneGeom (see VARARGIN)

% Automatic save defination after 3 minutes
t = timer;
t.StartDelay = 300;
t.Period = 300;
t.ExecutionMode = 'fixedRate';
handles.timer=t;
t.TimerFcn = {@SaveData_Callback, handles};

% Define the variables used to store geometries
data.crvs=[];  
data.chls=[]; 
data.isline=[]; 
data.hPanel=[];
data.n=0;  
data.new=0;  
data.dcrvs=[];  
data.disline=[];  
data.hLNum=[];
data.hPNum=[];
data.PlotLine='off';
data.PlotPoint='off';
data.filename=[]; 
set(handles.axes1, 'UserData', data); 
axis equal;

% Get mouse coordinates into CurrentCoords
set(gcf, 'WindowButtonDownFcn', {@SetCurrentCoords, handles, {'none'}});

% Choose default command line output for PlaneGeom
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PlaneGeom wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PlaneGeom_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function FilesMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FilesMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function NewFileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to NewFileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1, 'UserData');  
delete(data.chls);
data.crvs=[];  
data.chls=[]; 
data.isline=[]; 
data.n=[];
data.new=1;  
set(handles.axes1, 'UserData', data);  
RefreshPlot_Callback(hObject, eventdata, handles);
SaveAsMenu_Callback(hObject, eventdata, handles);


% --------------------------------------------------------------------
function OpenMenu_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in OpenButton.
function OpenButton_Callback(hObject, eventdata, handles)
% hObject    handle to OpenButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
OpenNurbs_Callback(hObject, eventdata, handles);


% --------------------------------------------------------------------
function OpenNurl_Callback(hObject, eventdata, handles)
% hObject    handle to OpenNurl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
{'**.mat;',...
 'MATLAB Files (*.mat)';
   '*.mat','MAT-files (*.mat)'; ...
   '*.*',  'All Files (*.*)'}, ...
   'Pick a file');
if filename
    opdata=load(filename, 'sdata');
    crvs=opdata.sdata.crvs;
    isline=opdata.sdata.isline;
    n=numel(isline);
    chls=zeros(1,n);
    for i=1:n
        if strcmp(crvs(i).form, 'POINTS')
            chls(i)=pntplot(crvs(i).coefs(1:2), 0.5, '.');
        else
            chls(i) = nrlcrvplot(crvs(i),50);
        end
    end
    axis equal;

    % Save file name
    data=get(handles.axes1, 'UserData');  
    data.filename=filename;

    % Set the axes limite
    XLim=get(handles.axes1,'XLim');
    YLim=get(handles.axes1,'YLim');
    set(handles.editXLim,'String', num2str(XLim));
    set(handles.editYLim,'String', num2str(YLim));

    % Select curves for further operation
    set(chls, 'ButtonDownFcn', {@Selected, handles, {'Selected'}});

    % Set Tag to line handls
    n=length(crvs);
    for i=1:n
        set(chls(i), 'Tag', num2str(i));
    end

    % Save data to 'UserData' of handles.axes1
    data.crvs=[data.crvs, crvs];  
    data.chls=[data.chls, chls]; 
    data.isline=[data.isline, isline]; 
    data.n=data.n+1;
    set(handles.axes1, 'UserData', data);  

    % Sart automatic saving
    n=data.n;
    if n==1
        start(handles.timer);
    end
end


% --- Executes on button press in OpenImage.
function OpenImage_Callback(hObject, eventdata, handles)
% hObject    handle to OpenImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, ~] = uigetfile( ...
{'*.jpg;*.tif;*.png;*.gif;*.bmp',...
 'Figure Files (*.jpg;*.tif;*.png;*.gif;*.bmp)';}, ...
   'Pick a file');
if filename
    A = imread(filename);
    h=image(A);

    % Define the variables used to store geometries
    data.crvs=[];  
    data.chls=[]; 
    data.isline=[]; 
    data.hPanel=[];
    data.n=0;  
    data.dcrvs=[];  
    data.disline=[];  
    data.hLNum=[];
    data.hPNum=[];
    data.PlotLine='off';
    data.PlotPoint='off';
    data.filename=[]; 
    data.image=h;
    set(handles.axes1, 'UserData', data); 
    
    % Get mouse coordinates into CurrentCoords
    set(gcf, 'WindowButtonDownFcn', {@SetCurrentCoords, handles, {'none'}});

    % Choose default command line output for PlaneGeom
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);
end


% --------------------------------------------------------------------
function SaveMenu_Callback(hObject, eventdata, handles)
% hObject    handle to SaveMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SaveData_Callback(hObject, eventdata, handles);


% --- Executes on button press in SaveData.
function SaveData_Callback(hObject, eventdata, handles)
% hObject    handle to SaveData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1,'UserData'); 
if isempty(data.filename) || data.new==1
    SaveAsMenu_Callback(hObject, eventdata, handles);
else
    sdata.crvs=data.crvs; 
    sdata.isline=data.isline; 
    save(data.filename, 'sdata');
end


% --------------------------------------------------------------------
function SaveAsMenu_Callback(hObject, eventdata, handles)
% hObject    handle to SaveAsMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uiputfile(...
 {'*.mat';'*.*'},...
 'Save as');
data=get(handles.axes1,'UserData'); 
sdata.crvs=data.crvs; 
sdata.isline=data.isline; 
if isequal(filename,0) || isequal(pathname,0)
else
    save(filename, 'sdata');
    data.filename=filename;
    set(handles.axes1,'UserData',data); 
end


% --------------------------------------------------------------------
function ExitMenu_Callback(hObject, eventdata, handles)
% hObject    handle to ExitMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stop(handles.timer);
button=questdlg('Are you sure to exit?','Exit','Yes','No','Yes');
if strcmp(button,'Yes')
    delete(handles.figure1);
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
stop(handles.timer);
button=questdlg('Are you sure to exit?','Exit','Yes','No','Yes');
if strcmp(button,'Yes')
    delete(hObject);
end


% --------------------------------------------------------------------
function CreateMenu_Callback(hObject, eventdata, handles)
% hObject    handle to CreateMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Points_Callback(hObject, eventdata, handles)
% hObject    handle to Points (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1, 'UserData');
if ~isempty(data.hPanel)
    delete(data.hPanel);
end
hp = uipanel('Title', 'Create Points', 'FontSize', 12, ...
             'BackgroundColor', 'white', ...
             'Position', [.78 .295 .2 .6], ...
             'DeleteFcn', {@DeleteFcn, handles, {'Mouse'}});
htext1 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Point:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 430 72 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit1=uicontrol('Parent', hp, 'Style', 'edit', 'String', '[0.0, 0.0]', ...
              'FontSize', 11, 'Position', [16 400 220 36], ...
              'HorizontalAlignment', 'left') ;
htext2 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'LineWidth:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 350 102 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit2=uicontrol('Parent', hp, 'Style', 'edit', 'String', '0.5', ...
              'FontSize', 11, 'Position', [16 320 220 36], ...
              'HorizontalAlignment', 'left') ;
htext3 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Marker:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 270 172 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit3=uicontrol('Parent', hp, 'Style', 'edit', 'String', '.', ...
              'FontSize', 11, 'Position', [16 240 220 36], ...
              'HorizontalAlignment', 'left') ;
argin{1}='Points'; 
argin{2}=hedit1; 
argin{3}=hedit2; 
argin{4}=hedit3; 
argin{5}=hp;
set(gcf, 'WindowButtonDownFcn', {@SetCurrentCoords, handles, argin});
set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, {'none'}}); 
set(data.chls, 'ButtonDownFcn', {@Selected, handles, {'none'}}); 
hPushCreate=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Create', 'FontSize', 11, ...
              'Position', [8 190 72 36], ...
              'Callback', {@Create, handles, argin}) ;
hPushReset=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Default', 'FontSize', 11, ...
              'Position', [98 190 72 36], ...
              'Callback', {@Reset, handles, argin}) ;
hPushCancel=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Cancel', 'FontSize', 11, ...
              'Position', [188 190 72 36], ...
              'Callback', {@Cancel, handles, argin}) ;
data.hPanel=hp;
set(handles.axes1, 'UserData', data); 


% --------------------------------------------------------------------
function CreateLines_Callback(hObject, eventdata, handles)
% hObject    handle to CreateLines (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Lines2Points_Callback(hObject, eventdata, handles)
% hObject    handle to Lines2Points (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1, 'UserData'); 
if ~isempty(data.hPanel)
    delete(data.hPanel);
end
hp = uipanel('Title', 'Create Lines', 'FontSize', 12, ...
             'BackgroundColor', 'white', ...
             'Position', [.78 .295 .2 .6], ...
             'DeleteFcn', {@DeleteFcn, handles, {'Mouse'}});
htext1 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Point 1:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 430 72 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit1=uicontrol('Parent', hp, 'Style', 'edit', 'String', '[0, 0]', ...
              'FontSize', 11, 'Position', [16 400 220 36], ...
              'HorizontalAlignment', 'left') ;
htext2 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Point 2:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 350 72 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit2=uicontrol('Parent', hp, 'Style', 'edit', 'String', '[1, 1]', ...
              'FontSize', 11, 'Position', [16 320 220 36], ...
              'HorizontalAlignment', 'left') ;
argin{1}='Line'; 
argin{2}=hedit1; 
argin{3}=hedit2; 
argin{4}=hp;
set(gcf, 'WindowButtonDownFcn', {@SetCurrentCoords, handles, argin});
set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, {'none'}}); 
set(data.chls, 'ButtonDownFcn', {@Selected, handles, {'none'}}); 
hPushCreate=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Create', 'FontSize', 11, ...
              'Position', [8 270 72 36], ...
              'Callback', {@Create, handles, argin}) ;
hPushReset=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Default', 'FontSize', 11, ...
              'Position', [98 270 72 36], ...
              'Callback', {@Reset, handles, argin}) ;
hPushCancel=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Cancel', 'FontSize', 11, ...
              'Position', [188 270 72 36], ...
              'Callback', {@Cancel, handles, argin}) ;
data.hPanel=hp;
set(handles.axes1, 'UserData', data); 


% --------------------------------------------------------------------
function Circles_Callback(hObject, eventdata, handles)
% hObject    handle to Circles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function CenterPoint_Callback(hObject, eventdata, handles)
% hObject    handle to CenterPoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1, 'UserData'); 
if ~isempty(data.hPanel)
    delete(data.hPanel);
end
hp = uipanel('Title', 'Create Circles', 'FontSize', 12, ...
             'BackgroundColor', 'white', ...
             'Position', [.78 .295 .2 .6], ...
             'DeleteFcn', {@DeleteFcn, handles, {'Mouse'}});
htext1 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Center:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 430 72 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit1=uicontrol('Parent', hp, 'Style', 'edit', 'String', '[0, 0]', ...
              'FontSize', 11, 'Position', [16 400 220 36], ...
              'HorizontalAlignment', 'left') ;
htext2 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Point:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 350 72 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit2=uicontrol('Parent', hp, 'Style', 'edit', 'String', '[1, 1]', ...
              'FontSize', 11, 'Position', [16 320 220 36], ...
              'HorizontalAlignment', 'left') ;
argin{1}='2PntsCircle'; 
argin{2}=hedit1; 
argin{3}=hedit2; 
argin{4}=hp;
set(gcf, 'WindowButtonDownFcn', {@SetCurrentCoords, handles, argin});
set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, {'none'}}); 
set(data.chls, 'ButtonDownFcn', {@Selected, handles, {'none'}}); 
hPushCreate=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Create', 'FontSize', 11, ...
              'Position', [8 270 72 36], ...
              'Callback', {@Create, handles, argin}) ;
hPushReset=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Default', 'FontSize', 11, ...
              'Position', [98 270 72 36], ...
              'Callback', {@Reset, handles, argin}) ;
hPushCancel=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Cancel', 'FontSize', 11, ...
              'Position', [188 270 72 36], ...
              'Callback', {@Cancel, handles, argin}) ;
data.hPanel=hp;
set(handles.axes1, 'UserData', data); 


% --------------------------------------------------------------------
function CenterRadius_Callback(hObject, eventdata, handles)
% hObject    handle to CenterPoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1, 'UserData'); 
if ~isempty(data.hPanel)
    delete(data.hPanel);
end
hp = uipanel('Title', 'Create Circles', 'FontSize', 12, ...
             'BackgroundColor', 'white', ...
             'Position', [.78 .295 .2 .6], ...
             'DeleteFcn', {@DeleteFcn, handles, {'Mouse'}});
htext1 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Center:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 430 72 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit1=uicontrol('Parent', hp, 'Style', 'edit', 'String', '[0, 0]', ...
              'FontSize', 11, 'Position', [16 400 220 36], ...
              'HorizontalAlignment', 'left') ;
htext2 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Radius:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 350 72 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit2=uicontrol('Parent', hp, 'Style', 'edit', 'String', '1.0', ...
              'FontSize', 11, 'Position', [16 320 220 36], ...
              'HorizontalAlignment', 'left') ;
argin{1}='CentRadius'; 
argin{2}=hedit1; 
argin{3}=hedit2; 
argin{4}=hp;
set(gcf, 'WindowButtonDownFcn', {@SetCurrentCoords, handles, argin});
set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, {'none'}}); 
set(data.chls, 'ButtonDownFcn', {@Selected, handles, {'none'}}); 
hPushCreate=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Create', 'FontSize', 11, ...
              'Position', [8 270 72 36], ...
              'Callback', {@Create, handles, argin}) ;
hPushReset=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Default', 'FontSize', 11, ...
              'Position', [98 270 72 36], ...
              'Callback', {@Reset, handles, argin}) ;
hPushCancel=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Cancel', 'FontSize', 11, ...
              'Position', [188 270 72 36], ...
              'Callback', {@Cancel, handles, argin}) ;
data.hPanel=hp;
set(handles.axes1, 'UserData', data); 


% --------------------------------------------------------------------
function Circles3pnts_Callback(hObject, eventdata, handles)
% hObject    handle to Circles3pnts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1, 'UserData'); 
if ~isempty(data.hPanel)
    delete(data.hPanel);
end
hp = uipanel('Title', 'Create Circles by 3 Points', 'FontSize', 12, ...
             'BackgroundColor', 'white', ...
             'Position', [.78 .295 .2 .6], ...
             'DeleteFcn', {@DeleteFcn, handles, {'Mouse'}});
htext1 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Point 1:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 430 72 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit1=uicontrol('Parent', hp, 'Style', 'edit', 'String', '[0.5, 0]', ...
              'FontSize', 11, 'Position', [16 400 220 36], ...
              'HorizontalAlignment', 'left') ;
htext2 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Point 2:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 350 72 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit2=uicontrol('Parent', hp, 'Style', 'edit', 'String', '[1, 1]', ...
              'FontSize', 11, 'Position', [16 320 220 36], ...
              'HorizontalAlignment', 'left') ;
htext3 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Point 3:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 270 172 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit3=uicontrol('Parent', hp, 'Style', 'edit', 'String', '[0, 1]', ...
              'FontSize', 11, 'Position', [16 240 220 36], ...
              'HorizontalAlignment', 'left') ;
argin{1}='Circles3pnts'; 
argin{2}=hedit1; 
argin{3}=hedit2; 
argin{4}=hedit3; 
argin{5}=hp;
set(gcf, 'WindowButtonDownFcn', {@SetCurrentCoords, handles, argin});
set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, {'none'}}); 
set(data.chls, 'ButtonDownFcn', {@Selected, handles, {'none'}}); 
hPushCreate=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Create', 'FontSize', 11, ...
              'Position', [8 190 72 36], ...
              'Callback', {@Create, handles, argin}) ;
hPushReset=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Default', 'FontSize', 11, ...
              'Position', [98 190 72 36], ...
              'Callback', {@Reset, handles, argin}) ;
hPushCancel=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Cancel', 'FontSize', 11, ...
              'Position', [188 190 72 36], ...
              'Callback', {@Cancel, handles, argin}) ;
data.hPanel=hp;
set(handles.axes1, 'UserData', data); 


% --------------------------------------------------------------------
function CtreateArcs_Callback(hObject, eventdata, handles)
% hObject    handle to CtreateArcs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function AcrCenterRadius_Callback(hObject, eventdata, handles)
% hObject    handle to AcrCenterRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1, 'UserData'); 
if ~isempty(data.hPanel)
    delete(data.hPanel);
end
hp = uipanel('Title', 'Create Arcs', 'FontSize', 12, ...
             'BackgroundColor', 'white', ...
             'Position', [.78 .295 .2 .6], ...
             'DeleteFcn', {@DeleteFcn, handles, {'Mouse'}});
htext1 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Radius:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 430 72 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit1=uicontrol('Parent', hp, 'Style', 'edit', 'String', '1', ...
              'FontSize', 11, 'Position', [16 400 220 36], ...
              'HorizontalAlignment', 'left') ;
htext2 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Center:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 350 72 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit2=uicontrol('Parent', hp, 'Style', 'edit', 'String', '[0, 0]', ...
              'FontSize', 11, 'Position', [16 320 220 36], ...
              'HorizontalAlignment', 'left') ;
htext3 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Start angle:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 270 172 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit3=uicontrol('Parent', hp, 'Style', 'edit', 'String', '0', ...
              'FontSize', 11, 'Position', [16 240 220 36], ...
              'HorizontalAlignment', 'left') ;
htext4 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'End angle:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 190 172 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit4=uicontrol('Parent', hp, 'Style', 'edit', 'String', '2*pi', ...
              'FontSize', 11, 'Position', [16 160 220 36], ...
              'HorizontalAlignment', 'left') ;
argin{1}='Circle'; 
argin{2}=hedit1; 
argin{3}=hedit2; 
argin{4}=hedit3; 
argin{5}=hedit4; 
argin{6}=hp;
set(gcf, 'WindowButtonDownFcn', {@SetCurrentCoords, handles, argin});
set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, {'none'}}); 
set(data.chls, 'ButtonDownFcn', {@Selected, handles, {'none'}}); 
hPushCreate=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Create', 'FontSize', 11, ...
              'Position', [8 110 72 36], ...
              'Callback', {@Create, handles, argin}) ;
hPushReset=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Default', 'FontSize', 11, ...
              'Position', [98 110 72 36], ...
              'Callback', {@Reset, handles, argin}) ;
hPushCancel=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Cancel', 'FontSize', 11, ...
              'Position', [188 110 72 36], ...
              'Callback', {@Cancel, handles, argin}) ;
data.hPanel=hp;
set(handles.axes1, 'UserData', data); 


% --------------------------------------------------------------------
function Arc3pnts_Callback(hObject, eventdata, handles)
% hObject    handle to Arc3pnts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1, 'UserData'); 
if ~isempty(data.hPanel)
    delete(data.hPanel);
end
hp = uipanel('Title', 'Create Arcs by 3 Points', 'FontSize', 12, ...
             'BackgroundColor', 'white', ...
             'Position', [.78 .295 .2 .6], ...
             'DeleteFcn', {@DeleteFcn, handles, {'Mouse'}});
htext1 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Point 1:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 430 72 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit1=uicontrol('Parent', hp, 'Style', 'edit', 'String', '[0.5, 0]', ...
              'FontSize', 11, 'Position', [16 400 220 36], ...
              'HorizontalAlignment', 'left') ;
htext2 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Point 2:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 350 72 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit2=uicontrol('Parent', hp, 'Style', 'edit', 'String', '[1, 1]', ...
              'FontSize', 11, 'Position', [16 320 220 36], ...
              'HorizontalAlignment', 'left') ;
htext3 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Point 3:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 270 172 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit3=uicontrol('Parent', hp, 'Style', 'edit', 'String', '[0, 1]', ...
              'FontSize', 11, 'Position', [16 240 220 36], ...
              'HorizontalAlignment', 'left') ;
argin{1}='Arc3pnts'; 
argin{2}=hedit1; 
argin{3}=hedit2; 
argin{4}=hedit3; 
argin{5}=hp;
set(gcf, 'WindowButtonDownFcn', {@SetCurrentCoords, handles, argin});
set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, {'none'}}); 
set(data.chls, 'ButtonDownFcn', {@Selected, handles, {'none'}}); 
hPushCreate=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Create', 'FontSize', 11, ...
              'Position', [8 190 72 36], ...
              'Callback', {@Create, handles, argin}) ;
hPushReset=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Default', 'FontSize', 11, ...
              'Position', [98 190 72 36], ...
              'Callback', {@Reset, handles, argin}) ;
hPushCancel=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Cancel', 'FontSize', 11, ...
              'Position', [188 190 72 36], ...
              'Callback', {@Cancel, handles, argin}) ;
data.hPanel=hp;
set(handles.axes1, 'UserData', data); 


% --------------------------------------------------------------------
function Triangles_Callback(hObject, eventdata, handles)
% hObject    handle to Triangles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1, 'UserData'); 
if ~isempty(data.hPanel)
    delete(data.hPanel);
end
hp = uipanel('Title', 'Create Triangles', 'FontSize', 12, ...
             'BackgroundColor', 'white', ...
             'Position', [.78 .295 .2 .6], ...
             'DeleteFcn', {@DeleteFcn, handles, {'Mouse'}});
htext1 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Point 1:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 430 72 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit1=uicontrol('Parent', hp, 'Style', 'edit', 'String', '[0.0, 0.0]', ...
              'FontSize', 11, 'Position', [16 400 220 36], ...
              'HorizontalAlignment', 'left') ;
htext2 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Point 2:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 350 72 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit2=uicontrol('Parent', hp, 'Style', 'edit', 'String', '[1.0, 0.0]', ...
              'FontSize', 11, 'Position', [16 320 220 36], ...
              'HorizontalAlignment', 'left') ;
htext3 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Point 3:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 270 172 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit3=uicontrol('Parent', hp, 'Style', 'edit', 'String', '[0.0, 1.0]', ...
              'FontSize', 11, 'Position', [16 240 220 36], ...
              'HorizontalAlignment', 'left') ;
argin{1}='Triangles'; 
argin{2}=hedit1; 
argin{3}=hedit2; 
argin{4}=hedit3; 
argin{5}=hp;
set(gcf, 'WindowButtonDownFcn', {@SetCurrentCoords, handles, argin});
set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, {'none'}}); 
set(data.chls, 'ButtonDownFcn', {@Selected, handles, {'none'}}); 
hPushCreate=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Create', 'FontSize', 11, ...
              'Position', [8 190 72 36], ...
              'Callback', {@Create, handles, argin}) ;
hPushReset=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Default', 'FontSize', 11, ...
              'Position', [98 190 72 36], ...
              'Callback', {@Reset, handles, argin}) ;
hPushCancel=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Cancel', 'FontSize', 11, ...
              'Position', [188 190 72 36], ...
              'Callback', {@Cancel, handles, argin}) ;
data.hPanel=hp;
set(handles.axes1, 'UserData', data); 


% --------------------------------------------------------------------
function Quadrangles_Callback(hObject, eventdata, handles)
% hObject    handle to Quadrangles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1, 'UserData'); 
if ~isempty(data.hPanel)
    delete(data.hPanel);
end
hp = uipanel('Title', 'Create Quadrangles', 'FontSize', 12, ...
             'BackgroundColor', 'white', ...
             'Position', [.78 .295 .2 .6], ...
             'DeleteFcn', {@DeleteFcn, handles, {'Mouse'}});
htext1 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Point 1:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 430 72 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit1=uicontrol('Parent', hp, 'Style', 'edit', 'String', '[0.0, 0.0]', ...
              'FontSize', 11, 'Position', [16 400 220 36], ...
              'HorizontalAlignment', 'left') ;
htext2 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Point 2:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 350 72 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit2=uicontrol('Parent', hp, 'Style', 'edit', 'String', '[1.0, 0.0]', ...
              'FontSize', 11, 'Position', [16 320 220 36], ...
              'HorizontalAlignment', 'left') ;
htext3 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Point 3:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 270 172 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit3=uicontrol('Parent', hp, 'Style', 'edit', 'String', '[0.0, 1.0]', ...
              'FontSize', 11, 'Position', [16 240 220 36], ...
              'HorizontalAlignment', 'left') ;
htext4 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Point 4:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 190 172 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit4=uicontrol('Parent', hp, 'Style', 'edit', 'String', '[1.0 1.0]', ...
              'FontSize', 11, 'Position', [16 160 220 36], ...
              'HorizontalAlignment', 'left') ;
argin{1}='Quadrangles'; 
argin{2}=hedit1; 
argin{3}=hedit2; 
argin{4}=hedit3; 
argin{5}=hedit4; 
argin{6}=hp;
set(gcf, 'WindowButtonDownFcn', {@SetCurrentCoords, handles, argin});
set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, {'none'}}); 
set(data.chls, 'ButtonDownFcn', {@Selected, handles, {'none'}}); 
hPushCreate=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Create', 'FontSize', 11, ...
              'Position', [8 110 72 36], ...
              'Callback', {@Create, handles, argin}) ;
hPushReset=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Default', 'FontSize', 11, ...
              'Position', [98 110 72 36], ...
              'Callback', {@Reset, handles, argin}) ;
hPushCancel=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Cancel', 'FontSize', 11, ...
              'Position', [188 110 72 36], ...
              'Callback', {@Cancel, handles, argin}) ;
data.hPanel=hp;
set(handles.axes1, 'UserData', data); 


% --------------------------------------------------------------------
function Polygons_Callback(hObject, eventdata, handles)
% hObject    handle to Polygons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1, 'UserData'); 
if ~isempty(data.hPanel)
    delete(data.hPanel);
end
hp = uipanel('Title', 'Create Polygons', 'FontSize', 12, ...
             'BackgroundColor', 'white', ...
             'Position', [.78 .295 .2 .6], ...
             'DeleteFcn', {@DeleteFcn, handles, {'Mouse'}});
htext1 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Number of Edges:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 430 172 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit1=uicontrol('Parent', hp, 'Style', 'edit', 'String', '5', ...
              'FontSize', 11, 'Position', [16 400 220 36], ...
              'HorizontalAlignment', 'left') ;
htext2 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Radius of Circumcircle:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 350 192 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit2=uicontrol('Parent', hp, 'Style', 'edit', 'String', '1.0', ...
              'FontSize', 11, 'Position', [16 320 220 36], ...
              'HorizontalAlignment', 'left') ;
htext3 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Center:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 270 172 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit3=uicontrol('Parent', hp, 'Style', 'edit', 'String', '[0.0, 0.0]', ...
              'FontSize', 11, 'Position', [16 240 220 36], ...
              'HorizontalAlignment', 'left') ;
htext4 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Angle of Lowest Edge:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 190 182 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit4=uicontrol('Parent', hp, 'Style', 'edit', 'String', '0.0', ...
              'FontSize', 11, 'Position', [16 160 220 36], ...
              'HorizontalAlignment', 'left') ;
argin{1}='Polygons'; 
argin{2}=hedit1; 
argin{3}=hedit2; 
argin{4}=hedit3; 
argin{5}=hedit4; 
argin{6}=hp;
set(gcf, 'WindowButtonDownFcn', {@SetCurrentCoords, handles, argin});
set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, {'none'}}); 
set(data.chls, 'ButtonDownFcn', {@Selected, handles, {'none'}}); 
hPushCreate=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Create', 'FontSize', 11, ...
              'Position', [8 110 72 36], ...
              'Callback', {@Create, handles, argin}) ;
hPushReset=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Default', 'FontSize', 11, ...
              'Position', [98 110 72 36], ...
              'Callback', {@Reset, handles, argin}) ;
hPushCancel=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Cancel', 'FontSize', 11, ...
              'Position', [188 110 72 36], ...
              'Callback', {@Cancel, handles, argin}) ;
data.hPanel=hp;
set(handles.axes1, 'UserData', data); 


% --------------------------------------------------------------------
function Analytical_Callback(hObject, eventdata, handles)
% hObject    handle to Analytical (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1, 'UserData'); 
if ~isempty(data.hPanel)
    delete(data.hPanel);
end
hp = uipanel('Title', 'Create Analytical Curves', 'FontSize', 12, ...
             'BackgroundColor', 'white', ...
             'Position', [.78 .295 .2 .6], ...
             'DeleteFcn', {@DeleteFcn, handles, {'Mouse'}});
htext1 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The order of B-Splines:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 430 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit1=uicontrol('Parent', hp, 'Style', 'edit', 'String', '2', ...
              'FontSize', 11, 'Position', [16 400 220 36], ...
              'HorizontalAlignment', 'left') ;
htext2 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Number of Control Points:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 350 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit2=uicontrol('Parent', hp, 'Style', 'edit', 'String', '15', ...
              'FontSize', 11, 'Position', [16 320 220 36], ...
              'HorizontalAlignment', 'left') ;
htext3 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Interval of Variable (s):', ...
              'BackgroundColor', 'white', ...
              'Position', [8 270 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit3=uicontrol('Parent', hp, 'Style', 'edit', 'String', '[0, 2*pi]', ...
              'FontSize', 11, 'Position', [16 240 220 36], ...
              'HorizontalAlignment', 'left') ;
htext4 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Function Handle of x:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 190 182 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit4=uicontrol('Parent', hp, 'Style', 'edit', 'String', '@(s)  4*cos(s)', ...
              'FontSize', 11, 'Position', [16 160 220 36], ...
              'HorizontalAlignment', 'left') ;
htext5 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Function Handle of y:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 110 182 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit5=uicontrol('Parent', hp, 'Style', 'edit', 'String', '@(s)  3*sin(s)', ...
              'FontSize', 11, 'Position', [16 80 220 36], ...
              'HorizontalAlignment', 'left') ;
argin{1}='Analytical'; 
argin{2}=hedit1; 
argin{3}=hedit2; 
argin{4}=hedit3; 
argin{5}=hedit4; 
argin{6}=hedit5; 
argin{7}=hp;
set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, {'none'}}); 
set(data.chls, 'ButtonDownFcn', {@Selected, handles, {'none'}}); 
hPushCreate=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Create', 'FontSize', 11, ...
              'Position', [8 30 72 36], ...
              'Callback', {@Create, handles, argin}) ;
hPushReset=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Default', 'FontSize', 11, ...
              'Position', [98 30 72 36], ...
              'Callback', {@Reset, handles, argin}) ;
hPushCancel=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Cancel', 'FontSize', 11, ...
              'Position', [188 30 72 36], ...
              'Callback', {@Cancel, handles, argin}) ;
data.hPanel=hp;
set(handles.axes1, 'UserData', data); 


% --------------------------------------------------------------------
function Splines_Callback(hObject, eventdata, handles)
% hObject    handle to Splines (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1, 'UserData'); 
if ~isempty(data.hPanel)
    delete(data.hPanel);
end
hp = uipanel('Title', 'Create Splines', 'FontSize', 12, ...
             'BackgroundColor', 'white', ...
             'Position', [.78 .295 .2 .6], ...
             'DeleteFcn', {@DeleteFcn, handles, {'Mouse'}});
htext1 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The order of B-Splines:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 430 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit1=uicontrol('Parent', hp, 'Style', 'edit', 'String', '2', ...
              'FontSize', 11, 'Position', [16 400 220 36], ...
              'HorizontalAlignment', 'left') ;
htext2 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The x-y Coordinates:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 350 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
cnames = {'x-Data','Y-Data'};
htable = uitable('Parent', hp, ...
              'ColumnWidth', {86}, ...
              'ColumnName', cnames, ...
              'ColumnEditable', [true true], ...
              'Position', [8 80 252 275], 'FontSize', 11);
argin{1}='Splines'; 
argin{2}=hedit1; 
argin{3}=htable; 
argin{4}=hp;
set(gcf, 'WindowButtonDownFcn', {@SetCurrentCoords, handles, argin});
set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, {'none'}}); 
set(data.chls, 'ButtonDownFcn', {@Selected, handles, {'none'}}); 
hPushCreate=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Create', 'FontSize', 11, ...
              'Position', [8 30 72 36], ...
              'Callback', {@Create, handles, argin}) ;
hPushReset=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Reset', 'FontSize', 11, ...
              'Position', [98 30 72 36], ...
              'Callback', {@Reset, handles, argin}) ;
hPushCancel=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Cancel', 'FontSize', 11, ...
              'Position', [188 30 72 36], ...
              'Callback', {@Cancel, handles, argin}) ;
data.hPanel=hp;
set(handles.axes1, 'UserData', data); 


function CurrentCoords_Callback(hObject, eventdata, handles)
% hObject    handle to CurrentCoords (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CurrentCoords as text
%        str2double(get(hObject,'String')) returns contents of CurrentCoords as a double


% --- Executes during object creation, after setting all properties.
function CurrentCoords_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CurrentCoords (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% Set the coordinates of InputCoords edit text
function SetCurrentCoords(hObject, callbackdata, handles, argin) 
    pos=get(handles.axes1,'Currentpoint');
    pos=pos(2,1:2);
    pos=num2str(pos);
    set(handles.CurrentCoords, 'String',pos); 
    switch argin{1}
        case 'Points'
            hedit1=argin{2}; 
            set(hedit1, 'String',pos); 
        case 'Line'
            hedit1=argin{2};
            hedit2=argin{3};
            [~, status1]=str2num(get(hedit1, 'String')); 
            [~, status2]=str2num(get(hedit2, 'String')); 
            if ~status1
                set(hedit1, 'String',pos); 
            elseif ~status2
                set(hedit2, 'String',pos); 
            end
        case '2PntsCircle'
            hedit1=argin{2}; 
            hedit2=argin{3}; 
            [~, status1]=str2num(get(hedit1, 'String')); 
            [~, status2]=str2num(get(hedit2, 'String')); 
            if ~status1
                set(hedit1, 'String',pos); 
            elseif ~status2
                set(hedit2, 'String',pos); 
            end
        case 'CentRadius'
            hedit1=argin{2}; 
            [~, status1]=str2num(get(hedit1, 'String')); 
            if ~status1
                set(hedit1, 'String',pos); 
            end
        case 'Circle'
            hedit2=argin{3}; 
            set(hedit2, 'String',pos); 
        case {'Circles3pnts', 'Arc3pnts', 'Triangles'}
            hedit1=argin{2}; 
            hedit2=argin{3}; 
            hedit3=argin{4}; 
            [~, status1]=str2num(get(hedit1, 'String')); 
            [~, status2]=str2num(get(hedit2, 'String')); 
            [~, status3]=str2num(get(hedit3, 'String')); 
            if ~status1
                set(hedit1, 'String',pos); 
            elseif ~status2
                set(hedit2, 'String',pos); 
            elseif ~status3
                set(hedit3, 'String',pos); 
            end
        case 'Quadrangles'
            hedit1=argin{2}; 
            hedit2=argin{3}; 
            hedit3=argin{4}; 
            hedit4=argin{5}; 
            [~, status1]=str2num(get(hedit1, 'String')); 
            [~, status2]=str2num(get(hedit2, 'String')); 
            [~, status3]=str2num(get(hedit3, 'String')); 
            [~, status4]=str2num(get(hedit4, 'String')); 
            if ~status1
                set(hedit1, 'String',pos); 
            elseif ~status2
                set(hedit2, 'String',pos); 
            elseif ~status3
                set(hedit3, 'String',pos); 
            elseif ~status4
                set(hedit4, 'String',pos); 
            end
        case 'Polygons'
            hedit3=argin{4}; 
            set(hedit3, 'String', pos); 
        case 'Splines'
            htable=argin{3};
            data=htable.Data;
            pos=str2num(pos);
            data=[data; pos];
            set(htable, 'Data', data);
        case {'LinePN', 'LinePT'}
            hedit1=argin{4}; 
            [~, status1]=str2num(get(hedit1, 'String')); 
            if ~status1
                set(hedit1, 'String',pos); 
            end
        case 'Rotate'
            hedit1=argin{5}; 
            [~, status1]=str2num(get(hedit1, 'String')); 
            if ~status1
                set(hedit1, 'String',pos); 
            end
    end


% Creat curves
function Create(hObject, callbackdata, handles, argin) 
    switch argin{1}
        case 'Points'
            hedit1=argin{2}; 
            hedit2=argin{3}; 
            hedit3=argin{4}; 
            pnt=str2num(get(hedit1, 'String')); 
            width=str2num(get(hedit2, 'String')); 
            maker=get(hedit3, 'String'); 
            crvs=pntmake(pnt);
            chls=pntplot(pnt, width, maker);
        case 'Line' 
            hedit1=argin{2}; 
            hedit2=argin{3}; 
            pnt1=str2num(get(hedit1, 'String')); 
            pnt2=str2num(get(hedit2, 'String')); 
            pnt1=GetCoords(handles, pnt1);
            pnt2=GetCoords(handles, pnt2);
            crvs=nrlline(pnt1, pnt2); 
            chls= nrlcrvplot(crvs,50); 
            set(hedit1, 'String', num2str(pnt2));            
            set(hedit2, 'String', num2str(rand(1,2)));
        case '2PntsCircle'
            hedit1=argin{2}; 
            hedit2=argin{3}; 
            pnt1=str2num(get(hedit1, 'String')); 
            pnt2=str2num(get(hedit2, 'String')); 
            pnt1=GetCoords(handles, pnt1);
            pnt2=GetCoords(handles, pnt2);
            R=sqrt(sum((pnt2-pnt1).^2));
            crvs=nrlcirc(R, [pnt1,0]);
            chls= nrlcrvplot(crvs,50); 
            set(hedit1, 'String', num2str(pnt1));            
            set(hedit2, 'String', num2str(pnt2));
        case 'CentRadius'
            hedit1=argin{2}; 
            hedit2=argin{3}; 
            pnt1=str2num(get(hedit1, 'String')); 
            R=str2num(get(hedit2, 'String')); 
            pnt1=GetCoords(handles, pnt1);
            crvs=nrlcirc(R, [pnt1,0]);
            chls= nrlcrvplot(crvs,50); 
            set(hedit1, 'String', num2str(pnt1)); 
        case 'Circle'
            hedit1=argin{2}; 
            hedit2=argin{3}; 
            hedit3=argin{4}; 
            hedit4=argin{5}; 
            radius=str2num(get(hedit1, 'String')); 
            center=str2num(get(hedit2, 'String')); 
            center=GetCoords(handles, center);
            sang=str2num(get(hedit3, 'String')); 
            eang=str2num(get(hedit4, 'String')); 
            crvs=nrlcirc(radius,center,sang,eang);
            chls = nrlcrvplot(crvs,50); 
        case 'Circles3pnts'
            hedit1=argin{2}; 
            hedit2=argin{3}; 
            hedit3=argin{4}; 
            pnt1=str2num(get(hedit1, 'String')); 
            pnt2=str2num(get(hedit2, 'String')); 
            pnt3=str2num(get(hedit3, 'String')); 
            pnt1=GetCoords(handles, pnt1);
            pnt2=GetCoords(handles, pnt2);
            pnt3=GetCoords(handles, pnt3); 
            crvs=nrl3ptscirc(pnt1, pnt2, pnt3);
            chls = nrlcrvplot(crvs,50);
        case 'Arc3pnts'
            hedit1=argin{2}; 
            hedit2=argin{3}; 
            hedit3=argin{4}; 
            pnt1=str2num(get(hedit1, 'String')); 
            pnt2=str2num(get(hedit2, 'String')); 
            pnt3=str2num(get(hedit3, 'String')); 
            pnt1=GetCoords(handles, pnt1);
            pnt2=GetCoords(handles, pnt2);
            pnt3=GetCoords(handles, pnt3); 
            crvs=nrl3pntsarc(pnt1, pnt2, pnt3);
            chls = nrlcrvplot(crvs,50);
        case 'Triangles'
            hedit1=argin{2}; 
            hedit2=argin{3}; 
            hedit3=argin{4}; 
            pnt1=str2num(get(hedit1, 'String')); 
            pnt2=str2num(get(hedit2, 'String')); 
            pnt3=str2num(get(hedit3, 'String')); 
            pnt1=GetCoords(handles, pnt1);
            pnt2=GetCoords(handles, pnt2);
            pnt3=GetCoords(handles, pnt3); 
            pnt1=[pnt1, 0];
            pnt2=[pnt2, 0];
            pnt3=[pnt3, 0];
            pnts={pnt1, pnt2, pnt3};
            crvs=nrltriangle(pnts);
            chls = nrlcrvplot(crvs,50);
        case 'Quadrangles'
            hedit1=argin{2}; 
            hedit2=argin{3}; 
            hedit3=argin{4}; 
            hedit4=argin{5}; 
            pnt1=str2num(get(hedit1, 'String')); 
            pnt2=str2num(get(hedit2, 'String')); 
            pnt3=str2num(get(hedit3, 'String')); 
            pnt4=str2num(get(hedit4, 'String')); 
            pnt1=GetCoords(handles, pnt1);
            pnt2=GetCoords(handles, pnt2);
            pnt3=GetCoords(handles, pnt3); 
            pnt4=GetCoords(handles, pnt4); 
            pnt1=[pnt1, 0];
            pnt2=[pnt2, 0];
            pnt3=[pnt3, 0];
            pnt4=[pnt4, 0];
            pnts={pnt1, pnt2, pnt3, pnt4};
            crvs=nrlquadr(pnts);
            chls = nrlcrvplot(crvs,50);
        case 'Polygons'
            hedit1=argin{2}; 
            hedit2=argin{3}; 
            hedit3=argin{4}; 
            hedit4=argin{5}; 
            n=str2num(get(hedit1, 'String')); 
            radius=str2num(get(hedit2, 'String')); 
            center=str2num(get(hedit3, 'String')); 
            center=GetCoords(handles, center); 
            ang=str2num(get(hedit4, 'String')); 
            crvs=nrlpolygon(n, radius, center, ang);
            chls = nrlcrvplot(crvs,50);
        case 'Analytical'
            hedit1=argin{2}; 
            hedit2=argin{3}; 
            hedit3=argin{4}; 
            hedit4=argin{5}; 
            hedit5=argin{6}; 
            p=str2num(get(hedit1, 'String')); 
            N=str2num(get(hedit2, 'String')); 
            a12=str2num(get(hedit3, 'String')); 
            funx=str2func(get(hedit4, 'String')); 
            funy=str2func(get(hedit5, 'String')); 
            a1=a12(1); a2=a12(2);
            crvs=nrlanalytical(p, N, a1, a2, funx, funy);
            chls = nrlcrvplot(crvs,50);
        case 'Splines'
            hedit1=argin{2}; 
            htable=argin{3}; 
            p=str2num(get(hedit1, 'String')); 
            data=htable.Data;
            x = data(:,1);
            y = data(:,2);
            crvs=nrlspline(p, x, y);
            chls = nrlcrvplot(crvs,50);
            Reset(hObject, callbackdata, handles, argin);
    end
    
    % Set the axes limite
    XLim=get(handles.axes1,'XLim');
    YLim=get(handles.axes1,'YLim');
    set(handles.editXLim,'String', num2str(XLim));
    set(handles.editYLim,'String', num2str(YLim));
    
    % Select curves for further operation
    set(chls, 'ButtonDownFcn', {@Selected, handles, {'Selected'}});
    
    % Judge whether the entity is a line or a point
    data=get(handles.axes1, 'UserData');  
    n=length(chls);
    if strcmp(argin{1},  'Points')
        isline=0;
    else
        isline=ones(1,n);
    end
    
    % Set Tag to line handls
    m=length(data.crvs);
    for i=1:n
        set(chls(i), 'Tag', num2str(m+i));
    end
    
    % Save data to 'UserData' of handles.axes1
    data.crvs=[data.crvs, crvs];  
    data.chls=[data.chls, chls]; 
    data.isline=[data.isline, isline]; 
    data.n=data.n+1;
    set(handles.axes1, 'UserData', data);  
    
    % Sart automatic saving
    n=data.n;
    if n==1
        start(handles.timer);
    end

% Reset data
function Reset(hObject, callbackdata, handles, argin) 
    switch argin{1}
        case 'Points'
            hedit1=argin{2}; 
            hedit2=argin{3}; 
            hedit3=argin{4}; 
            set(hedit1, 'String', '[0.0 0.0]');
            set(hedit2, 'String', '0.5');
            set(hedit3, 'String', '.');
        case {'Line', '2PntsCircle'}
            hedit1=argin{2}; 
            hedit2=argin{3}; 
            set(hedit1, 'String', '[0, 0]');
            set(hedit2, 'String', '[1, 1]');
        case 'CentRadius'
            hedit1=argin{2}; 
            hedit2=argin{3}; 
            set(hedit1, 'String', '[0, 0]');
            set(hedit2, 'String', '1.0');
        case 'Circle'
            hedit1=argin{2}; 
            hedit2=argin{3}; 
            hedit3=argin{4}; 
            hedit4=argin{5}; 
            set(hedit1, 'String', '1');
            set(hedit2, 'String', '[0, 0]');
            set(hedit3, 'String', '0');
            set(hedit4, 'String', '2*pi');
        case {'Circles3pnts', 'Arc3pnts'}
            hedit1=argin{2}; 
            hedit2=argin{3}; 
            hedit3=argin{4}; 
            set(hedit1, 'String', '[0.5, 0]');
            set(hedit2, 'String', '[1, 1]');
            set(hedit3, 'String', '[0, 1]');
        case 'Triangles'
            hedit1=argin{2}; 
            hedit2=argin{3}; 
            hedit3=argin{4}; 
            set(hedit1, 'String', '[0.0, 0.0]');
            set(hedit2, 'String', '[1.0, 0.0]');
            set(hedit3, 'String', '[0.0, 1.0]');
        case 'Quadrangles'
            hedit1=argin{2}; 
            hedit2=argin{3}; 
            hedit3=argin{4}; 
            hedit4=argin{5}; 
            set(hedit1, 'String', '[0.0, 0.0]');
            set(hedit2, 'String', '[1.0, 0.0]');
            set(hedit3, 'String', '[0.0, 1.0]');
            set(hedit4, 'String', '[1.0, 1.0]');
        case 'Polygons'
            hedit1=argin{2}; 
            hedit2=argin{3}; 
            hedit3=argin{4}; 
            hedit4=argin{5}; 
            set(hedit1, 'String', '5');
            set(hedit2, 'String', '1.0');
            set(hedit3, 'String', '[0.0, 0.0]');
            set(hedit4, 'String', '0.0');
        case 'Analytical'
            hedit1=argin{2}; 
            hedit2=argin{3}; 
            hedit3=argin{4}; 
            hedit4=argin{5}; 
            hedit5=argin{6}; 
            set(hedit1, 'String', '2');
            set(hedit2, 'String', '15');
            set(hedit3, 'String', '[0, 2*pi]');
            set(hedit4, 'String', '@(s)  4*cos(s)');
            set(hedit5, 'String', '@(s)  3*sin(s)');
        case 'Splines'
            htable=argin{3};
            set(htable, 'Data', []);
        case {'Deleted', 'Extract', 'Center'}
            htable=argin{2};
            set(htable, 'Data', []);
        case 'Save'
            htable=argin{2};
            hedit1=argin{4};
            set(htable, 'Data', []);
            set(hedit1, 'String', 'savedata');
        case 'Fillet'
            htable=argin{2};
            hedit=argin{4}; 
            set(htable, 'Data', []);
            set(hedit, 'String', '0.25');
        case 'Measure'
            htable=argin{2};
            hedit=argin{4}; 
            set(htable, 'Data', []);
            set(hedit, 'String', '...');
        case 'FilletTA'
            htable=argin{2};
            hedit=argin{4}; 
            set(htable, 'Data', []);
            set(hedit, 'String', '[1 1]');
        case {'Cut', 'Mirror'}
            htable=argin{2};
            hedit=argin{4};            
            set(htable, 'Data', []);
            set(hedit,  'String', 'input a curve number');
        case 'Glue'
            htable=argin{2};
            set(htable, 'Data', []);
        case 'Split'
            htable=argin{2};
            hedit=argin{4};            
            set(htable, 'Data', []);
            set(hedit,  'String', '[0.25, 0.5, 0.75]');
        case 'Translate'
            htable=argin{2};
            hedit=argin{4};      
            hedit2=argin{5};   
            set(htable, 'Data', []);
            set(hedit,  'String', '[1.0, 0.5]');
            set(hedit2,  'String', '1');
        case 'Scale'
             htable=argin{2};
            hedit=argin{4};      
            hedit2=argin{5};   
            set(htable, 'Data', []);
            set(hedit,  'String', '[0.5, 0.5]');
            set(hedit2,  'String', '1');
        case 'Double'
            htable=argin{2};
            set(htable, 'Data', []);
        case 'Rotate'
            htable=argin{2};
            hedit=argin{4};      
            hedit2=argin{5};   
            hedit3=argin{6};   
            set(htable, 'Data', []);
            set(hedit,  'String', 'pi/2');
            set(hedit2,  'String', '[0.0, 0.0]');
            set(hedit3,  'String', '1');
        case 'CutByPoint'
            htable=argin{2};
            hedit=argin{4};            
            set(htable, 'Data', []);
            set(hedit,  'String', 'input a point number');
            
            data=get(handles.axes1,'UserData'); 
            set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, argin}); 
        case {'LinePN', 'LinePT'}
            hedit2=argin{2};
            hedit1=argin{4};        
            hedit3=argin{5};       
            set(hedit1,  'String', '[0, 0]');
            set(hedit2,  'String', 'input a line number');
            set(hedit3,  'String', '1.0');         
        case 'RePoints'
            hedit1=argin{2};            
            set(hedit1,  'String', '0');
    end


% Cancle creating createmenu
function Cancel(hObject, callbackdata, handles, argin) 
    switch argin{1}
        case 'Points'
            hp=argin{5};            
        case {'Line', 'Splines', '2PntsCircle', 'CentRadius'}
            hp=argin{4};            
        case {'Circle', 'Quadrangles', 'Polygons'}
            hp=argin{6};
        case {'Circles3pnts', 'Arc3pnts', 'Triangles'}
            hp=argin{5};
        case 'Analytical'
            hp=argin{7};
        case {'Deleted', 'Fillet', 'FilletTA', 'Cut', 'Extract', 'Center', 'LinePN', 'LinePT', 'RePoints'}            
            hp=argin{3};
        case {'Measure', 'CutByPoint', 'Split', 'Glue', 'Translate', 'Rotate', 'Mirror', 'Double', 'Save'}
            hp=argin{3};
        case 'Scale'
            hp=argin{3};
    end
    delete(hp);
    set(gcf, 'WindowButtonDownFcn', {@SetCurrentCoords, handles, {'none'}});


% Make a point of this package
function crvs=pntmake(pnt)
crvs.form='POINTS';
crvs.dim=3;
crvs.number=1;
crvs.coefs=[pnt, 0];
crvs.knots=1;
crvs.intervals=[0,1];
crvs.order=1;


% Plot a point of this package
function chls=pntplot(pnt, width, maker) 
chls = line(pnt(1), pnt(2), 'LineWidth', width, 'Marker', maker);


% Get a point from the data storing all curves and points
function pnt=GetCoords(handles, pnt) 
n=length(pnt);
if n==1
    data=get(handles.axes1,'UserData');
    num=data.PNumber(pnt);    
    if ~data.isline(num)
        crvs=data.crvs(num);
        pnt=crvs.coefs(1:2);
    end
end


function editXLim_Callback(hObject, eventdata, handles)
% hObject    handle to editXLim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editXLim as text
%        str2double(get(hObject,'String')) returns contents of editXLim as a double


% --- Executes during object creation, after setting all properties.
function editXLim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editXLim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editYLim_Callback(hObject, eventdata, handles)
% hObject    handle to editYLim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editYLim as text
%        str2double(get(hObject,'String')) returns contents of editYLim as a double


% --- Executes during object creation, after setting all properties.
function editYLim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editYLim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SetXLim.
function SetXLim_Callback(hObject, eventdata, handles)
% hObject    handle to SetXLim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
XLim=get(handles.editXLim,'String');
[XLim, status]=str2num(XLim);
if status
    xlim(XLim);
end


% --- Executes on button press in SetYLim.
function SetYLim_Callback(hObject, eventdata, handles)
% hObject    handle to SetYLim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
YLim=get(handles.editYLim,'String');
[YLim, status]=str2num(YLim);
if status
    ylim(YLim);
end


% --- Executes on button press in AutoAxisButton.
function AutoAxisButton_Callback(hObject, eventdata, handles)
% hObject    handle to AutoAxisButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axis auto; 
axis equal;


% --------------------------------------------------------------------
function PlotMenu_Callback(hObject, eventdata, handles)
% hObject    handle to PlotMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Generate_Callback(hObject, eventdata, handles)
% hObject    handle to Generate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1,'UserData');
isline=data.isline; 
n=length(isline);
nL=length(find(isline==1));
nP=n-nL;
s=1; t=1;
LNumber=zeros(1, nL);
PNumber=zeros(1, nP);
T2LNum=zeros(1,n);
T2PNum=zeros(1,n);
for i=1:n
    if isline(i)
        LNumber(s)=i;
        s=s+1;
    else
        PNumber(t)=i;
        t=t+1;
    end
    set(data.chls(i), 'Tag', num2str(i));
end
for i=1:nL
    T2LNum(LNumber(i))=i;
end
for i=1:nP
    T2PNum(PNumber(i))=i;
end
data.LNumber=LNumber;
data.PNumber=PNumber;
data.T2LNum=T2LNum;
data.T2PNum=T2PNum;
set(handles.axes1, 'UserData', data);


% --------------------------------------------------------------------
function LineNumber_Callback(hObject, eventdata, handles)
% hObject    handle to LineNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Remove the old line numbers
data=get(handles.axes1,'UserData'); 
if ~isempty(data.hLNum)
    m=length(data.hLNum);
    for i=1:m
        if ishandle(data.hLNum(i))
            delete(data.hLNum(i));
        end
    end
end

% Generate and plot new line numbers
Generate_Callback(hObject, eventdata, handles); 
data=get(handles.axes1,'UserData'); 
data.PlotLine='on';
crvs=data.crvs; 
LNumber=data.LNumber; 
n=length(LNumber);
hLNum=zeros(1,n);
for i=1:n
    p = nrleval(crvs(LNumber(i)), 0.5);
    hLNum(i) = text(p(1), p(2), num2str(i), 'Color', 'r');
end
data.hLNum=hLNum;
set(handles.axes1, 'UserData', data);

% --------------------------------------------------------------------
function PointNumber_Callback(hObject, eventdata, handles)
% hObject    handle to PointNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Remove the old point numbers
data=get(handles.axes1,'UserData'); 
if ~isempty(data.hPNum)
    m=length(data.hPNum);
    for i=1:m
        if ishandle(data.hPNum(i))
            delete(data.hPNum(i));
        end
    end
end

% Generate and plot new point numbers
Generate_Callback(hObject, eventdata, handles); 
data=get(handles.axes1,'UserData'); 
data.PlotPoint='on';
PNumber=data.PNumber; 
crvs=data.crvs; 
n=length(PNumber);
hPNum=zeros(1,n);
for i=1:n
    p = crvs(PNumber(i)).coefs;
    hPNum(i) = text(p(1), p(2), num2str(i), 'Color', 'b');
    set(hPNum(i), 'Tag', num2str(i));
end

% Select points for further operation
set(hPNum, 'ButtonDownFcn', {@Selected, handles, {'Selected'}});

% Save data
data.hPNum=hPNum;
set(handles.axes1, 'UserData', data);


% --- Executes on button press in RefreshPlot.
function RefreshPlot_Callback(hObject, eventdata, handles)
% hObject    handle to RefreshPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Generate_Callback(hObject, eventdata, handles); 

data=get(handles.axes1,'UserData'); 
if strcmp(data.PlotLine, 'off')
    if ~isempty(data.hLNum)
        m=length(data.hLNum);
        for i=1:m
            if ishandle(data.hLNum(i))
                delete(data.hLNum(i));
            end
        end
        data.hLNum=[];
        set(handles.axes1, 'UserData', data);
    end
elseif strcmp(data.PlotLine, 'on')
    LineNumber_Callback(hObject, eventdata, handles);
end

if strcmp(data.PlotPoint, 'off')
    if ~isempty(data.hPNum)
        m=length(data.hPNum);
        for i=1:m
            if ishandle(data.hPNum(i))
                delete(data.hPNum(i));
            end
        end
        data.hPNum=[];
        set(handles.axes1, 'UserData', data);
    end
elseif strcmp(data.PlotPoint, 'on')
    PointNumber_Callback(hObject, eventdata, handles);
end


% --------------------------------------------------------------------
function LineNumberOff_Callback(hObject, eventdata, handles)
% hObject    handle to LineNumberOff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1,'UserData'); 
data.PlotLine='off';
if ~isempty(data.hLNum)
    m=length(data.hLNum);
    for i=1:m
        if ishandle(data.hLNum(i))
            delete(data.hLNum(i));
        end
    end
    data.hLNum=[];
end
set(handles.axes1, 'UserData', data);


% --------------------------------------------------------------------
function PointNumberOff_Callback(hObject, eventdata, handles)
% hObject    handle to PointNumberOff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1,'UserData'); 
data.PlotPoint='off';
if ~isempty(data.hPNum)
    m=length(data.hPNum);
    for i=1:m
        if ishandle(data.hPNum(i))
            delete(data.hPNum(i));
        end
    end
    data.hPNum=[];
end
set(handles.axes1, 'UserData', data);


% --------------------------------------------------------------------
function EditMenu_Callback(hObject, eventdata, handles)
% hObject    handle to EditMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function DeleteCurves_Callback(hObject, eventdata, handles)
% hObject    handle to DeleteCurves (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1, 'UserData'); 
if ~isempty(data.hPanel)
    delete(data.hPanel);
end
Generate_Callback(hObject, eventdata, handles);
hp = uipanel('Title', 'Delete Curves and Points', 'FontSize', 12, ...
             'BackgroundColor', 'white', ...
             'Position', [.78 .295 .2 .6], ...
             'DeleteFcn', {@DeleteFcn, handles, {'Selected'}});
htext = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The curves to be deleted:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 430 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
htable = uitable('Parent', hp, ...
              'ColumnWidth', {86}, ...
              'ColumnName', {'Types', 'Number'}, ...
              'ColumnEditable', [false true], ...
              'Position', [8 80 252 340], 'FontSize', 11);
argin{1}='Deleted'; 
argin{2}=htable; 
argin{3}=hp;
data=get(handles.axes1,'UserData'); 
set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, argin});
set(data.chls, 'ButtonDownFcn', {@Selected, handles, argin});
hPushCreate=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Confirm', 'FontSize', 11, ...
              'Position', [8 30 72 36], ...
              'Callback', {@Edit, handles, argin}) ;
hPushReset=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Reset', 'FontSize', 11, ...
              'Position', [98 30 72 36], ...
              'Callback', {@Reset, handles, argin}) ;
hPushCancel=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Cancel', 'FontSize', 11, ...
              'Position', [188 30 72 36], ...
              'Callback', {@Cancel, handles, argin}) ;
data.hPanel=hp;
set(handles.axes1, 'UserData', data); 


% --------------------------------------------------------------------
function SaveLines_Callback(hObject, eventdata, handles)
% hObject    handle to SaveLines (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1, 'UserData'); 
if ~isempty(data.hPanel)
    delete(data.hPanel);
end
Generate_Callback(hObject, eventdata, handles);
hp = uipanel('Title', 'Save Curves and Points', 'FontSize', 12, ...
             'BackgroundColor', 'white', ...
             'Position', [.78 .295 .2 .6], ...
             'DeleteFcn', {@DeleteFcn, handles, {'Selected'}});
htext1 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The file name:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 430 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit1=uicontrol('Parent', hp, 'Style', 'edit', 'String', 'savedata', ...
              'FontSize', 11, 'Position', [16 400 220 36], ...
              'HorizontalAlignment', 'left') ;    
htext = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The curves to be saved:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 350 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');      
htable = uitable('Parent', hp, ...
              'ColumnWidth', {86}, ...
              'ColumnName', {'Types', 'Number'}, ...
              'ColumnEditable', [false true], ...
              'Position', [8 80 252 275], 'FontSize', 11);
argin{1}='Save'; 
argin{2}=htable; 
argin{3}=hp;
argin{4}=hedit1;
data=get(handles.axes1,'UserData'); 
set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, argin});
set(data.chls, 'ButtonDownFcn', {@Selected, handles, argin});
hPushCreate=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Confirm', 'FontSize', 11, ...
              'Position', [8 30 72 36], ...
              'Callback', {@Edit, handles, argin}) ;
hPushReset=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Reset', 'FontSize', 11, ...
              'Position', [98 30 72 36], ...
              'Callback', {@Reset, handles, argin}) ;
hPushCancel=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Cancel', 'FontSize', 11, ...
              'Position', [188 30 72 36], ...
              'Callback', {@Cancel, handles, argin}) ;
data.hPanel=hp;
set(handles.axes1, 'UserData', data); 


% --------------------------------------------------------------------
function FilletLines_Callback(hObject, eventdata, handles)
% hObject    handle to FilletLines (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1, 'UserData'); 
if ~isempty(data.hPanel)
    delete(data.hPanel);
end
Generate_Callback(hObject, eventdata, handles);
hp = uipanel('Title', 'Fillet two lines', 'FontSize', 12, ...
             'BackgroundColor', 'white', ...
             'Position', [.78 .295 .2 .6], ...
             'DeleteFcn', {@DeleteFcn, handles, {'Selected'}});
htext1 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The radius of the fillet:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 430 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit1=uicontrol('Parent', hp, 'Style', 'edit', 'String', '0.25', ...
              'FontSize', 11, 'Position', [16 400 220 36], ...
              'HorizontalAlignment', 'left') ;
htext = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The two lines to be filleted:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 350 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
htable = uitable('Parent', hp, ...
              'ColumnWidth', {86}, ...
              'ColumnName', {'Types', 'Number'}, ...
              'ColumnEditable', [false true], ...
              'Position', [8 80 252 260], 'FontSize', 11);
argin{1}='Fillet'; 
argin{2}=htable; 
argin{3}=hp; 
argin{4}=hedit1; 
data=get(handles.axes1,'UserData'); 
set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, argin}); 
set(data.chls, 'ButtonDownFcn', {@Selected, handles, argin}); 
hPushCreate=uicontrol('Parent', hp, 'Style', 'pushbutton', ... 
              'String', 'Confirm', 'FontSize', 11, ... 
              'Position', [8 30 72 36], ... 
              'Callback', {@Edit, handles, argin}) ; 
hPushReset=uicontrol('Parent', hp, 'Style', 'pushbutton', ... 
              'String', 'Reset', 'FontSize', 11, ... 
              'Position', [98 30 72 36], ...
              'Callback', {@Reset, handles, argin}) ;
hPushCancel=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Cancel', 'FontSize', 11, ...
              'Position', [188 30 72 36], ...
              'Callback', {@Cancel, handles, argin}) ;
data.hPanel=hp;
set(handles.axes1, 'UserData', data); 


% --------------------------------------------------------------------
function Measure_Callback(hObject, eventdata, handles)
% hObject    handle to Measure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1, 'UserData'); 
if ~isempty(data.hPanel)
    delete(data.hPanel);
end
Generate_Callback(hObject, eventdata, handles);
hp = uipanel('Title', 'Measure a curve', 'FontSize', 12, ...
             'BackgroundColor', 'white', ...
             'Position', [.78 .295 .2 .6], ...
             'DeleteFcn', {@DeleteFcn, handles, {'Selected'}});
htext1 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The length of the curve:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 430 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit1=uicontrol('Parent', hp, 'Style', 'edit', 'String', '...', ...
              'FontSize', 11, 'Position', [16 400 220 36], ...
              'HorizontalAlignment', 'left') ;
htext = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The curve to be measured:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 350 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
htable = uitable('Parent', hp, ...
              'ColumnWidth', {86}, ...
              'ColumnName', {'Types', 'Number'}, ...
              'ColumnEditable', [false true], ...
              'Position', [8 80 252 260], 'FontSize', 11);
argin{1}='Measure'; 
argin{2}=htable; 
argin{3}=hp; 
argin{4}=hedit1; 
data=get(handles.axes1,'UserData'); 
set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, argin}); 
set(data.chls, 'ButtonDownFcn', {@Selected, handles, argin}); 
hPushCreate=uicontrol('Parent', hp, 'Style', 'pushbutton', ... 
              'String', 'Confirm', 'FontSize', 11, ... 
              'Position', [8 30 72 36], ... 
              'Callback', {@Edit, handles, argin}) ; 
hPushReset=uicontrol('Parent', hp, 'Style', 'pushbutton', ... 
              'String', 'Reset', 'FontSize', 11, ... 
              'Position', [98 30 72 36], ...
              'Callback', {@Reset, handles, argin}) ;
hPushCancel=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Cancel', 'FontSize', 11, ...
              'Position', [188 30 72 36], ...
              'Callback', {@Cancel, handles, argin}) ;
data.hPanel=hp;
set(handles.axes1, 'UserData', data); 


% --------------------------------------------------------------------
function TangentArcs_Callback(hObject, eventdata, handles)
% hObject    handle to TangentArcs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1, 'UserData'); 
if ~isempty(data.hPanel)
    delete(data.hPanel);
end
Generate_Callback(hObject, eventdata, handles);
hp = uipanel('Title', 'Tagent arcs', 'FontSize', 12, ...
             'BackgroundColor', 'white', ...
             'Position', [.78 .295 .2 .6], ...
             'DeleteFcn', {@DeleteFcn, handles, {'Selected'}});
htext1 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Small or large one?', ...
              'BackgroundColor', 'white', ...
              'Position', [8 430 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit1=uicontrol('Parent', hp, 'Style', 'edit', 'String', '[1 1]', ...
              'FontSize', 11, 'Position', [16 400 220 36], ...
              'HorizontalAlignment', 'left') ;
htext = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The two lines to be filleted:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 350 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
htable = uitable('Parent', hp, ...
              'ColumnWidth', {86}, ...
              'ColumnName', {'Types', 'Number'}, ...
              'ColumnEditable', [false true], ...
              'Position', [8 80 252 260], 'FontSize', 11);
argin{1}='FilletTA'; 
argin{2}=htable; 
argin{3}=hp; 
argin{4}=hedit1; 
data=get(handles.axes1,'UserData'); 
set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, argin}); 
set(data.chls, 'ButtonDownFcn', {@Selected, handles, argin}); 
hPushCreate=uicontrol('Parent', hp, 'Style', 'pushbutton', ... 
              'String', 'Confirm', 'FontSize', 11, ... 
              'Position', [8 30 72 36], ... 
              'Callback', {@Edit, handles, argin}) ; 
hPushReset=uicontrol('Parent', hp, 'Style', 'pushbutton', ... 
              'String', 'Reset', 'FontSize', 11, ... 
              'Position', [98 30 72 36], ...
              'Callback', {@Reset, handles, argin}) ;
hPushCancel=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Cancel', 'FontSize', 11, ...
              'Position', [188 30 72 36], ...
              'Callback', {@Cancel, handles, argin}) ;
data.hPanel=hp;
set(handles.axes1, 'UserData', data); 


% --------------------------------------------------------------------
function CutCurvesMenu_Callback(hObject, eventdata, handles)
% hObject    handle to CutCurvesMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function CutByCurve_Callback(hObject, eventdata, handles)
% hObject    handle to CutByCurve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1, 'UserData'); 
if ~isempty(data.hPanel)
    delete(data.hPanel);
end
Generate_Callback(hObject, eventdata, handles);
hp = uipanel('Title', 'Cut curves', 'FontSize', 12, ...
             'BackgroundColor', 'white', ...
             'Position', [.78 .295 .2 .6], ...
             'DeleteFcn', {@DeleteFcn, handles, {'Selected'}});
htext1 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The curve to cut:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 430 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit1=uicontrol('Parent', hp, 'Style', 'edit', 'String', 'input a curve number', ...
              'FontSize', 11, 'Position', [16 400 220 36], ...
              'HorizontalAlignment', 'left') ;
htext = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The curves to be cutted:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 350 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
htable = uitable('Parent', hp, ...
              'ColumnWidth', {86}, ...
              'ColumnName', {'Types', 'Number'}, ...
              'ColumnEditable', [false true], ...
              'Position', [8 80 252 275], 'FontSize', 11);
argin{1}='Cut'; 
argin{2}=htable; 
argin{3}=hp; 
argin{4}=hedit1;
data=get(handles.axes1,'UserData'); 
set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, argin}); 
set(data.chls, 'ButtonDownFcn', {@Selected, handles, argin}); 
hPushCreate=uicontrol('Parent', hp, 'Style', 'pushbutton', ... 
              'String', 'Confirm', 'FontSize', 11, ... 
              'Position', [8 30 72 36], ... 
              'Callback', {@Edit, handles, argin}) ; 
hPushReset=uicontrol('Parent', hp, 'Style', 'pushbutton', ... 
              'String', 'Reset', 'FontSize', 11, ... 
              'Position', [98 30 72 36], ...
              'Callback', {@Reset, handles, argin}) ;
hPushCancel=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Cancel', 'FontSize', 11, ...
              'Position', [188 30 72 36], ...
              'Callback', {@Cancel, handles, argin}) ;
data.hPanel=hp;
set(handles.axes1, 'UserData', data); 


% --------------------------------------------------------------------
function CutByPoint_Callback(hObject, eventdata, handles)
% hObject    handle to CutByPoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1, 'UserData'); 
if ~isempty(data.hPanel)
    delete(data.hPanel);
end
Generate_Callback(hObject, eventdata, handles);
hp = uipanel('Title', 'Cut Curve by a Point', 'FontSize', 12, ...
             'BackgroundColor', 'white', ...
             'Position', [.78 .295 .2 .6], ...
             'DeleteFcn', {@DeleteFcn, handles, {'Selected'}});
htext1 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The point to cut:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 430 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit1=uicontrol('Parent', hp, 'Style', 'edit', 'String', 'input a point number', ...
              'FontSize', 11, 'Position', [16 400 220 36], ...
              'HorizontalAlignment', 'left') ;
htext = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The curve to be cutted:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 350 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
htable = uitable('Parent', hp, ...
              'ColumnWidth', {86}, ...
              'ColumnName', {'Types', 'Number'}, ...
              'ColumnEditable', [false true], ...
              'Position', [8 80 252 275], 'FontSize', 11);
argin{1}='CutByPoint'; 
argin{2}=htable; 
argin{3}=hp; 
argin{4}=hedit1;
data=get(handles.axes1,'UserData'); 
set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, argin}); 
set(data.chls, 'ButtonDownFcn', {@Selected, handles, argin}); 
hPushCreate=uicontrol('Parent', hp, 'Style', 'pushbutton', ... 
              'String', 'Confirm', 'FontSize', 11, ... 
              'Position', [8 30 72 36], ... 
              'Callback', {@Edit, handles, argin}) ; 
hPushReset=uicontrol('Parent', hp, 'Style', 'pushbutton', ... 
              'String', 'Reset', 'FontSize', 11, ... 
              'Position', [98 30 72 36], ...
              'Callback', {@Reset, handles, argin}) ;
hPushCancel=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Cancel', 'FontSize', 11, ...
              'Position', [188 30 72 36], ...
              'Callback', {@Cancel, handles, argin}) ;
data.hPanel=hp;
set(handles.axes1, 'UserData', data); 


% --------------------------------------------------------------------
function SplitCurves_Callback(hObject, eventdata, handles)
% hObject    handle to SplitCurves (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1, 'UserData'); 
if ~isempty(data.hPanel)
    delete(data.hPanel);
end
Generate_Callback(hObject, eventdata, handles);
hp = uipanel('Title', 'Split curves', 'FontSize', 12, ...
             'BackgroundColor', 'white', ...
             'Position', [.78 .295 .2 .6], ...
             'DeleteFcn', {@DeleteFcn, handles, {'Selected'}});
htext1 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Knots to split (or number):', ...
              'BackgroundColor', 'white', ...
              'Position', [8 430 222 56], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit1=uicontrol('Parent', hp, 'Style', 'edit', 'String', '[0.25, 0.5, 0.75]', ...
              'FontSize', 11, 'Position', [16 400 220 36], ...
              'HorizontalAlignment', 'left') ;
htext = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The curves to be splitted:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 350 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
htable = uitable('Parent', hp, ...
              'ColumnWidth', {86}, ...
              'ColumnName', {'Types', 'Number'}, ...
              'ColumnEditable', [false true], ...
              'Position', [8 80 252 275], 'FontSize', 11);
argin{1}='Split'; 
argin{2}=htable; 
argin{3}=hp; 
argin{4}=hedit1;
data=get(handles.axes1,'UserData'); 
set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, argin}); 
set(data.chls, 'ButtonDownFcn', {@Selected, handles, argin}); 
hPushCreate=uicontrol('Parent', hp, 'Style', 'pushbutton', ... 
              'String', 'Confirm', 'FontSize', 11, ... 
              'Position', [8 30 72 36], ... 
              'Callback', {@Edit, handles, argin}) ; 
hPushReset=uicontrol('Parent', hp, 'Style', 'pushbutton', ... 
              'String', 'Reset', 'FontSize', 11, ... 
              'Position', [98 30 72 36], ...
              'Callback', {@Reset, handles, argin}) ;
hPushCancel=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Cancel', 'FontSize', 11, ...
              'Position', [188 30 72 36], ...
              'Callback', {@Cancel, handles, argin}) ;
data.hPanel=hp;
set(handles.axes1, 'UserData', data); 


% --------------------------------------------------------------------
function GlueCurves_Callback(hObject, eventdata, handles)
% hObject    handle to GlueCurves (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1, 'UserData'); 
if ~isempty(data.hPanel)
    delete(data.hPanel);
end
Generate_Callback(hObject, eventdata, handles);
hp = uipanel('Title', 'Glue curves', 'FontSize', 12, ...
             'BackgroundColor', 'white', ...
             'Position', [.78 .295 .2 .6], ...
             'DeleteFcn', {@DeleteFcn, handles, {'Selected'}});
htext = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The curves to be glued:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 430 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
htable = uitable('Parent', hp, ...
              'ColumnWidth', {86}, ...
              'ColumnName', {'Types', 'Number'}, ...
              'ColumnEditable', [false true], ...
              'Position', [8 80 252 355], 'FontSize', 11);
argin{1}='Glue'; 
argin{2}=htable; 
argin{3}=hp; 
data=get(handles.axes1,'UserData'); 
set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, argin}); 
set(data.chls, 'ButtonDownFcn', {@Selected, handles, argin}); 
hPushCreate=uicontrol('Parent', hp, 'Style', 'pushbutton', ... 
              'String', 'Confirm', 'FontSize', 11, ... 
              'Position', [8 30 72 36], ... 
              'Callback', {@Edit, handles, argin}) ; 
hPushReset=uicontrol('Parent', hp, 'Style', 'pushbutton', ... 
              'String', 'Reset', 'FontSize', 11, ... 
              'Position', [98 30 72 36], ...
              'Callback', {@Reset, handles, argin}) ;
hPushCancel=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Cancel', 'FontSize', 11, ...
              'Position', [188 30 72 36], ...
              'Callback', {@Cancel, handles, argin}) ;
data.hPanel=hp;
set(handles.axes1, 'UserData', data); 


% --------------------------------------------------------------------
function ExtractMenu_Callback(hObject, eventdata, handles)
% hObject    handle to ExtractMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ExtractPoints_Callback(hObject, eventdata, handles)
% hObject    handle to ExtractPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1, 'UserData'); 
if ~isempty(data.hPanel)
    delete(data.hPanel);
end
Generate_Callback(hObject, eventdata, handles);
hp = uipanel('Title', 'Extract points', 'FontSize', 12, ...
             'BackgroundColor', 'white', ...
             'Position', [.78 .295 .2 .6], ...
             'DeleteFcn', {@DeleteFcn, handles, {'Selected'}});
htext = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The curves to be extracted:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 430 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
htable = uitable('Parent', hp, ...
              'ColumnWidth', {86}, ...
              'ColumnName', {'Types', 'Number'}, ...
              'ColumnEditable', [false true], ...
              'Position', [8 80 252 340], 'FontSize', 11);
argin{1}='Extract'; 
argin{2}=htable; 
argin{3}=hp; 
data=get(handles.axes1,'UserData'); 
set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, argin}); 
set(data.chls, 'ButtonDownFcn', {@Selected, handles, argin}); 
hPushCreate=uicontrol('Parent', hp, 'Style', 'pushbutton', ... 
              'String', 'Confirm', 'FontSize', 11, ... 
              'Position', [8 30 72 36], ... 
              'Callback', {@Edit, handles, argin}) ; 
hPushReset=uicontrol('Parent', hp, 'Style', 'pushbutton', ... 
              'String', 'Reset', 'FontSize', 11, ... 
              'Position', [98 30 72 36], ...
              'Callback', {@Reset, handles, argin}) ;
hPushCancel=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Cancel', 'FontSize', 11, ...
              'Position', [188 30 72 36], ...
              'Callback', {@Cancel, handles, argin}) ;
data.hPanel=hp;
set(handles.axes1, 'UserData', data); 


% --------------------------------------------------------------------
function ExtractCenters_Callback(hObject, eventdata, handles)
% hObject    handle to ExtractCenters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1, 'UserData'); 
if ~isempty(data.hPanel)
    delete(data.hPanel);
end
Generate_Callback(hObject, eventdata, handles);
hp = uipanel('Title', 'Extract points', 'FontSize', 12, ...
             'BackgroundColor', 'white', ...
             'Position', [.78 .295 .2 .6], ...
             'DeleteFcn', {@DeleteFcn, handles, {'Selected'}});
htext = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The curves to be extracted:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 430 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
htable = uitable('Parent', hp, ...
              'ColumnWidth', {86}, ...
              'ColumnName', {'Types', 'Number'}, ...
              'ColumnEditable', [false true], ...
              'Position', [8 80 252 340], 'FontSize', 11);
argin{1}='Center'; 
argin{2}=htable; 
argin{3}=hp; 
data=get(handles.axes1,'UserData'); 
set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, argin}); 
set(data.chls, 'ButtonDownFcn', {@Selected, handles, argin}); 
hPushCreate=uicontrol('Parent', hp, 'Style', 'pushbutton', ... 
              'String', 'Confirm', 'FontSize', 11, ... 
              'Position', [8 30 72 36], ... 
              'Callback', {@Edit, handles, argin}) ; 
hPushReset=uicontrol('Parent', hp, 'Style', 'pushbutton', ... 
              'String', 'Reset', 'FontSize', 11, ... 
              'Position', [98 30 72 36], ...
              'Callback', {@Reset, handles, argin}) ;
hPushCancel=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Cancel', 'FontSize', 11, ...
              'Position', [188 30 72 36], ...
              'Callback', {@Cancel, handles, argin}) ;
data.hPanel=hp;
set(handles.axes1, 'UserData', data); 


% --------------------------------------------------------------------
function LinesPointNorm_Callback(hObject, eventdata, handles)
% hObject    handle to LinesPointNorm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1, 'UserData'); 
if ~isempty(data.hPanel)
    delete(data.hPanel);
end
Generate_Callback(hObject, eventdata, handles);
hp = uipanel('Title', 'Create Lines', 'FontSize', 12, ...
             'BackgroundColor', 'white', ...
             'Position', [.78 .295 .2 .6], ...
             'DeleteFcn', {@DeleteFcn, handles, {'Selected'}});
htext1 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The point:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 430 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit1=uicontrol('Parent', hp, 'Style', 'edit', 'String', '[0, 0]', ...
              'FontSize', 11, 'Position', [16 400 220 36], ...
              'HorizontalAlignment', 'left') ;
htext2 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The normal line:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 350 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit2=uicontrol('Parent', hp, 'Style', 'edit', 'String', 'input a line number', ...
              'FontSize', 11, 'Position', [16 320 220 36], ...
              'HorizontalAlignment', 'left') ;
htext3 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Length of the line:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 270 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit3=uicontrol('Parent', hp, 'Style', 'edit', 'String', '1.0', ...
              'FontSize', 11, 'Position', [16 240 220 36], ...
              'HorizontalAlignment', 'left') ;
argin{1}='LinePT'; 
argin{2}=hedit2; 
argin{3}=hp; 
argin{4}=hedit1; 
argin{5}=hedit3; 
data=get(handles.axes1,'UserData'); 
set(gcf, 'WindowButtonDownFcn', {@SetCurrentCoords, handles, argin});
set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, argin}); 
set(data.chls, 'ButtonDownFcn', {@Selected, handles, argin}); 
hPushCreate=uicontrol('Parent', hp, 'Style', 'pushbutton', ... 
              'String', 'Confirm', 'FontSize', 11, ... 
              'Position', [8 30 72 36], ... 
              'Callback', {@Edit, handles, argin}) ; 
hPushReset=uicontrol('Parent', hp, 'Style', 'pushbutton', ... 
              'String', 'Reset', 'FontSize', 11, ... 
              'Position', [98 30 72 36], ...
              'Callback', {@Reset, handles, argin}) ;
hPushCancel=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Cancel', 'FontSize', 11, ...
              'Position', [188 30 72 36], ...
              'Callback', {@Cancel, handles, argin}) ;
data.hPanel=hp;
set(handles.axes1, 'UserData', data); 


% --------------------------------------------------------------------
function LinesPointTang_Callback(hObject, eventdata, handles)
% hObject    handle to LinesPointTang (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1, 'UserData'); 
if ~isempty(data.hPanel)
    delete(data.hPanel);
end
Generate_Callback(hObject, eventdata, handles);
hp = uipanel('Title', 'Create Lines', 'FontSize', 12, ...
             'BackgroundColor', 'white', ...
             'Position', [.78 .295 .2 .6], ...
             'DeleteFcn', {@DeleteFcn, handles, {'Selected'}});
htext1 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The point:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 430 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit1=uicontrol('Parent', hp, 'Style', 'edit', 'String', '[0, 0]', ...
              'FontSize', 11, 'Position', [16 400 220 36], ...
              'HorizontalAlignment', 'left') ;
htext2 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The tangent line:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 350 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit2=uicontrol('Parent', hp, 'Style', 'edit', 'String', 'input a line number', ...
              'FontSize', 11, 'Position', [16 320 220 36], ...
              'HorizontalAlignment', 'left') ;
htext3 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Length of the line:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 270 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit3=uicontrol('Parent', hp, 'Style', 'edit', 'String', '1.0', ...
              'FontSize', 11, 'Position', [16 240 220 36], ...
              'HorizontalAlignment', 'left') ;
argin{1}='LinePN'; 
argin{2}=hedit2; 
argin{3}=hp; 
argin{4}=hedit1; 
argin{5}=hedit3; 
data=get(handles.axes1,'UserData'); 
set(gcf, 'WindowButtonDownFcn', {@SetCurrentCoords, handles, argin});
set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, argin}); 
set(data.chls, 'ButtonDownFcn', {@Selected, handles, argin}); 
hPushCreate=uicontrol('Parent', hp, 'Style', 'pushbutton', ... 
              'String', 'Confirm', 'FontSize', 11, ... 
              'Position', [8 30 72 36], ... 
              'Callback', {@Edit, handles, argin}) ; 
hPushReset=uicontrol('Parent', hp, 'Style', 'pushbutton', ... 
              'String', 'Reset', 'FontSize', 11, ... 
              'Position', [98 30 72 36], ...
              'Callback', {@Reset, handles, argin}) ;
hPushCancel=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Cancel', 'FontSize', 11, ...
              'Position', [188 30 72 36], ...
              'Callback', {@Cancel, handles, argin}) ;
data.hPanel=hp;
set(handles.axes1, 'UserData', data); 


% --------------------------------------------------------------------
function RemovePoints_Callback(hObject, eventdata, handles)
% hObject    handle to RemovePoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 data=get(handles.axes1, 'UserData'); 
if ~isempty(data.hPanel)
    delete(data.hPanel);
end
Generate_Callback(hObject, eventdata, handles);
hp = uipanel('Title', 'Remove Points', 'FontSize', 12, ...
             'BackgroundColor', 'white', ...
             'Position', [.78 .295 .2 .6], ...
             'DeleteFcn', {@DeleteFcn, handles, {'Selected'}});
htext1 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Only dumplicated?', ...
              'BackgroundColor', 'white', ...
              'Position', [8 430 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit1=uicontrol('Parent', hp, 'Style', 'edit', 'String', '0', ...
              'FontSize', 11, 'Position', [16 400 220 36], ...
              'HorizontalAlignment', 'left') ;
argin{1}='RePoints'; 
argin{2}=hedit1; 
argin{3}=hp;
data=get(handles.axes1,'UserData'); 
set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, argin});
set(data.chls, 'ButtonDownFcn', {@Selected, handles, argin});
hPushCreate=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Confirm', 'FontSize', 11, ...
              'Position', [8 30 72 36], ...
              'Callback', {@Edit, handles, argin}) ;
hPushReset=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Reset', 'FontSize', 11, ...
              'Position', [98 30 72 36], ...
              'Callback', {@Reset, handles, argin}) ;
hPushCancel=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Cancel', 'FontSize', 11, ...
              'Position', [188 30 72 36], ...
              'Callback', {@Cancel, handles, argin}) ;
data.hPanel=hp;
set(handles.axes1, 'UserData', data); 


% --------------------------------------------------------------------
function TransformMenu_Callback(hObject, eventdata, handles)
% hObject    handle to TransformMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Translation_Callback(hObject, eventdata, handles)
% hObject    handle to Translation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1, 'UserData'); 
if ~isempty(data.hPanel)
    delete(data.hPanel);
end
Generate_Callback(hObject, eventdata, handles);
hp = uipanel('Title', 'Translate curves', 'FontSize', 12, ...
             'BackgroundColor', 'white', ...
             'Position', [.78 .295 .2 .6], ...
             'DeleteFcn', {@DeleteFcn, handles, {'Selected'}});
htext1 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The translation vector:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 460 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit1=uicontrol('Parent', hp, 'Style', 'edit', 'String', '[1.0, 0.5]', ...
              'FontSize', 11, 'Position', [16 430 220 36], ...
              'HorizontalAlignment', 'left') ;
htext2 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Retain the original curves?', ...
              'BackgroundColor', 'white', ...
              'Position', [8 380 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit2=uicontrol('Parent', hp, 'Style', 'edit', 'String', '1', ...
              'FontSize', 11, 'Position', [16 350 220 36], ...
              'HorizontalAlignment', 'left') ;
htext = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The curves to be translated:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 300 282 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
htable = uitable('Parent', hp, ...
              'ColumnWidth', {86}, ...
              'ColumnName', {'Types', 'Number'}, ...
              'ColumnEditable', [false true], ...
              'Position', [8 80 252 225], 'FontSize', 11);
argin{1}='Translate'; 
argin{2}=htable; 
argin{3}=hp; 
argin{4}=hedit1;
argin{5}=hedit2;
data=get(handles.axes1,'UserData'); 
set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, argin}); 
set(data.chls, 'ButtonDownFcn', {@Selected, handles, argin}); 
hPushCreate=uicontrol('Parent', hp, 'Style', 'pushbutton', ... 
              'String', 'Confirm', 'FontSize', 11, ... 
              'Position', [8 30 72 36], ... 
              'Callback', {@Edit, handles, argin}) ; 
hPushReset=uicontrol('Parent', hp, 'Style', 'pushbutton', ... 
              'String', 'Reset', 'FontSize', 11, ... 
              'Position', [98 30 72 36], ...
              'Callback', {@Reset, handles, argin}) ;
hPushCancel=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Cancel', 'FontSize', 11, ...
              'Position', [188 30 72 36], ...
              'Callback', {@Cancel, handles, argin}) ;
data.hPanel=hp;
set(handles.axes1, 'UserData', data); 


% --------------------------------------------------------------------
function Rotate_Callback(hObject, eventdata, handles)
% hObject    handle to Rotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1, 'UserData'); 
if ~isempty(data.hPanel)
    delete(data.hPanel);
end
Generate_Callback(hObject, eventdata, handles);
hp = uipanel('Title', 'Rotate curves', 'FontSize', 12, ...
             'BackgroundColor', 'white', ...
             'Position', [.78 .295 .2 .6], ...
             'DeleteFcn', {@DeleteFcn, handles, {'Selected'}});
htext1 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The angle to rotate the curve:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 460 282 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit1=uicontrol('Parent', hp, 'Style', 'edit', 'String', 'pi/2', ...
              'FontSize', 11, 'Position', [16 430 220 36], ...
              'HorizontalAlignment', 'left') ;
htext2 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The reference point:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 380 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit2=uicontrol('Parent', hp, 'Style', 'edit', 'String', '[0.0, 0.0]', ...
              'FontSize', 11, 'Position', [16 350 220 36], ...
              'HorizontalAlignment', 'left') ;
htext3 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Retain the original curves?', ...
              'BackgroundColor', 'white', ...
              'Position', [8 300 272 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit3=uicontrol('Parent', hp, 'Style', 'edit', 'String', '1', ...
              'FontSize', 11, 'Position', [16 270 220 36], ...
              'HorizontalAlignment', 'left') ;
htext = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The curves to be rotated:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 220 300 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
htable = uitable('Parent', hp, ...
              'ColumnWidth', {86}, ...
              'ColumnName', {'Types', 'Number'}, ...
              'ColumnEditable', [false true], ...
              'Position', [8 80 252 145], 'FontSize', 11);
argin{1}='Rotate'; 
argin{2}=htable; 
argin{3}=hp; 
argin{4}=hedit1;
argin{5}=hedit2;
argin{6}=hedit3;
data=get(handles.axes1,'UserData'); 
set(gcf, 'WindowButtonDownFcn', {@SetCurrentCoords, handles, argin});
set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, argin}); 
set(data.chls, 'ButtonDownFcn', {@Selected, handles, argin}); 
hPushCreate=uicontrol('Parent', hp, 'Style', 'pushbutton', ... 
              'String', 'Confirm', 'FontSize', 11, ... 
              'Position', [8 30 72 36], ... 
              'Callback', {@Edit, handles, argin}) ; 
hPushReset=uicontrol('Parent', hp, 'Style', 'pushbutton', ... 
              'String', 'Reset', 'FontSize', 11, ... 
              'Position', [98 30 72 36], ...
              'Callback', {@Reset, handles, argin}) ;
hPushCancel=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Cancel', 'FontSize', 11, ...
              'Position', [188 30 72 36], ...
              'Callback', {@Cancel, handles, argin}) ;
data.hPanel=hp;
set(handles.axes1, 'UserData', data); 


% --------------------------------------------------------------------
function ScaleTrans_Callback(hObject, eventdata, handles)
% hObject    handle to ScaleTrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1, 'UserData'); 
if ~isempty(data.hPanel)
    delete(data.hPanel);
end
Generate_Callback(hObject, eventdata, handles);
hp = uipanel('Title', 'Scale curves', 'FontSize', 12, ...
             'BackgroundColor', 'white', ...
             'Position', [.78 .295 .2 .6], ...
             'DeleteFcn', {@DeleteFcn, handles, {'Selected'}});
htext1 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The scale:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 460 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit1=uicontrol('Parent', hp, 'Style', 'edit', 'String', '[0.5, 0.5]', ...
              'FontSize', 11, 'Position', [16 430 220 36], ...
              'HorizontalAlignment', 'left') ;
htext2 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Retain the original curves?', ...
              'BackgroundColor', 'white', ...
              'Position', [8 380 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit2=uicontrol('Parent', hp, 'Style', 'edit', 'String', '1', ...
              'FontSize', 11, 'Position', [16 350 220 36], ...
              'HorizontalAlignment', 'left') ;
htext = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The curves to be scaled:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 300 282 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
htable = uitable('Parent', hp, ...
              'ColumnWidth', {86}, ...
              'ColumnName', {'Types', 'Number'}, ...
              'ColumnEditable', [false true], ...
              'Position', [8 80 252 225], 'FontSize', 11);
argin{1}='Scale'; 
argin{2}=htable; 
argin{3}=hp; 
argin{4}=hedit1;
argin{5}=hedit2;
data=get(handles.axes1,'UserData'); 
set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, argin}); 
set(data.chls, 'ButtonDownFcn', {@Selected, handles, argin}); 
hPushCreate=uicontrol('Parent', hp, 'Style', 'pushbutton', ... 
              'String', 'Confirm', 'FontSize', 11, ... 
              'Position', [8 30 72 36], ... 
              'Callback', {@Edit, handles, argin}) ; 
hPushReset=uicontrol('Parent', hp, 'Style', 'pushbutton', ... 
              'String', 'Reset', 'FontSize', 11, ... 
              'Position', [98 30 72 36], ...
              'Callback', {@Reset, handles, argin}) ;
hPushCancel=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Cancel', 'FontSize', 11, ...
              'Position', [188 30 72 36], ...
              'Callback', {@Cancel, handles, argin}) ;
data.hPanel=hp;
set(handles.axes1, 'UserData', data); 


% --------------------------------------------------------------------
function Mirror_Callback(hObject, eventdata, handles)
% hObject    handle to Mirror (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1, 'UserData'); 
if ~isempty(data.hPanel)
    delete(data.hPanel);
end
Generate_Callback(hObject, eventdata, handles);
hp = uipanel('Title', 'Mirror curves', 'FontSize', 12, ...
             'BackgroundColor', 'white', ...
             'Position', [.78 .295 .2 .6], ...
             'DeleteFcn', {@DeleteFcn, handles, {'Selected'}});
htext1 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The reference line:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 430 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit1=uicontrol('Parent', hp, 'Style', 'edit', 'String', 'input a curve number', ...
              'FontSize', 11, 'Position', [16 400 220 36], ...
              'HorizontalAlignment', 'left') ;
htext = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The curves to be mirrored:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 350 222 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
htable = uitable('Parent', hp, ...
              'ColumnWidth', {86}, ...
              'ColumnName', {'Types', 'Number'}, ...
              'ColumnEditable', [false true], ...
              'Position', [8 80 252 275], 'FontSize', 11);
argin{1}='Mirror'; 
argin{2}=htable; 
argin{3}=hp; 
argin{4}=hedit1;
data=get(handles.axes1,'UserData'); 
set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, argin}); 
set(data.chls, 'ButtonDownFcn', {@Selected, handles, argin}); 
hPushCreate=uicontrol('Parent', hp, 'Style', 'pushbutton', ... 
              'String', 'Confirm', 'FontSize', 11, ... 
              'Position', [8 30 72 36], ... 
              'Callback', {@Edit, handles, argin}) ; 
hPushReset=uicontrol('Parent', hp, 'Style', 'pushbutton', ... 
              'String', 'Reset', 'FontSize', 11, ... 
              'Position', [98 30 72 36], ...
              'Callback', {@Reset, handles, argin}) ;
hPushCancel=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Cancel', 'FontSize', 11, ...
              'Position', [188 30 72 36], ...
              'Callback', {@Cancel, handles, argin}) ;
data.hPanel=hp;
set(handles.axes1, 'UserData', data); 


% --------------------------------------------------------------------
function Double_Callback(hObject, eventdata, handles)
% hObject    handle to Double (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1, 'UserData'); 
if ~isempty(data.hPanel)
    delete(data.hPanel);
end
Generate_Callback(hObject, eventdata, handles);
hp = uipanel('Title', 'Double curves', 'FontSize', 12, ...
             'BackgroundColor', 'white', ...
             'Position', [.78 .295 .2 .6], ...
             'DeleteFcn', {@DeleteFcn, handles, {'Selected'}});
htext1 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'Knots insertion:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 460 272 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit1=uicontrol('Parent', hp, 'Style', 'edit', 'String', '0', ...
              'FontSize', 11, 'Position', [16 430 220 36], ...
              'HorizontalAlignment', 'left') ;
htext2 = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The distance between curves:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 380 272 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
hedit2=uicontrol('Parent', hp, 'Style', 'edit', 'String', '0.1', ...
              'FontSize', 11, 'Position', [16 350 220 36], ...
              'HorizontalAlignment', 'left') ;
htext = uicontrol('Parent', hp, 'Style', 'text', 'String', 'The curves to be doubled:', ...
              'BackgroundColor', 'white', ...
              'Position', [8 300 282 36], 'FontSize', 11, ...
              'HorizontalAlignment', 'left');
htable = uitable('Parent', hp, ...
              'ColumnWidth', {86}, ...
              'ColumnName', {'Types', 'Number'}, ...
              'ColumnEditable', [false true], ...
              'Position', [8 80 252 225], 'FontSize', 11);
argin{1}='Double'; 
argin{2}=htable; 
argin{3}=hp; 
argin{4}=hedit1;
argin{5}=hedit2;
data=get(handles.axes1,'UserData'); 
set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, argin}); 
set(data.chls, 'ButtonDownFcn', {@Selected, handles, argin}); 
hPushCreate=uicontrol('Parent', hp, 'Style', 'pushbutton', ... 
              'String', 'Confirm', 'FontSize', 11, ... 
              'Position', [8 30 72 36], ... 
              'Callback', {@Edit, handles, argin}) ; 
hPushReset=uicontrol('Parent', hp, 'Style', 'pushbutton', ... 
              'String', 'Reset', 'FontSize', 11, ... 
              'Position', [98 30 72 36], ...
              'Callback', {@Reset, handles, argin}) ;
hPushCancel=uicontrol('Parent', hp, 'Style', 'pushbutton', ...
              'String', 'Cancel', 'FontSize', 11, ...
              'Position', [188 30 72 36], ...
              'Callback', {@Cancel, handles, argin}) ;
data.hPanel=hp;
set(handles.axes1, 'UserData', data); 


% Select curves for further operation
function Selected(hObject, callbackdata, handles, argin) 
    switch argin{1}
        case {'Selected', 'RePoints'}
            Type=get(callbackdata.Source, 'Type');
            cnum=str2num(get(callbackdata.Source, 'Tag'));
            data=get(handles.axes1,'UserData');
            if strcmp(Type, 'text')
                cnum=data.PNumber(cnum);
            end
            set(data.chls(cnum), 'LineWidth', 1.5, 'Color', 'm');
            pause(0.5);
            set(data.chls(cnum), 'LineWidth', 0.5, 'Color', 'k');
        case 'Deleted'
            Type=get(callbackdata.Source, 'Type');
            cnum=str2num(get(callbackdata.Source, 'Tag'));
            data=get(handles.axes1,'UserData');
            htable=argin{2};
            tdata=htable.Data;
            if strcmp(Type, 'text')
                pos={'Point', cnum};
                cnum=data.PNumber(cnum); 
            else
                pos={'Curve', data.T2LNum(cnum)};
            end
            tdata=[tdata; pos];
            tdata=removetabel(tdata);
            set(htable, 'Data', tdata); 
            set(data.chls(cnum), 'LineWidth', 1.5, 'Color', 'm');
            pause(0.5);
            set(data.chls(cnum), 'LineWidth', 0.5, 'Color', 'k');
        case {'Fillet', 'FilletTA'}
            Type=get(callbackdata.Source, 'Type');
            cnum=str2num(get(callbackdata.Source, 'Tag'));
            data=get(handles.axes1,'UserData');
            htable=argin{2};
            tdata=htable.Data;
            if strcmp(Type, 'line')
                pos={'Curve', data.T2LNum(cnum)};
                n=numel(tdata)/2;            
                if n<2
                    tdata=[tdata; pos];
                    tdata=removetabel(tdata);
                    set(htable, 'Data', tdata); 
                end
                set(data.chls(cnum), 'LineWidth', 1.5, 'Color', 'm');
                pause(0.5); 
                set(data.chls(cnum), 'LineWidth', 0.5, 'Color', 'k');                
            end
        case 'Measure'
            Type=get(callbackdata.Source, 'Type');
            cnum=str2num(get(callbackdata.Source, 'Tag'));
            data=get(handles.axes1,'UserData');
            htable=argin{2};
            tdata=htable.Data;
            if strcmp(Type, 'line')
                pos={'Curve', data.T2LNum(cnum)};
                n=numel(tdata)/2;            
                if n<1
                    tdata=[tdata; pos];
                    set(htable, 'Data', tdata); 
                end
                set(data.chls(cnum), 'LineWidth', 1.5, 'Color', 'm');
                pause(0.5); 
                set(data.chls(cnum), 'LineWidth', 0.5, 'Color', 'k');                
            end
        case 'Cut'
            htable=argin{2};
            hedit=argin{4};
            
            empt=str2num(get(hedit, 'String')); 
            cnum=str2num(get(callbackdata.Source, 'Tag'));
            data=get(handles.axes1,'UserData');            
            tdata=htable.Data;
            if strcmp(data.crvs(cnum).form, 'L-NURL')
                if isempty(empt)
                    lnum=data.T2LNum(cnum);
                    set(hedit, 'String', num2str(lnum));
                    set(data.chls(cnum), 'LineWidth', 1.5, 'Color', 'm');
                    pause(0.5); 
                    set(data.chls(cnum), 'LineWidth', 0.5, 'Color', 'k'); 
                else
                    pos={'Curve', data.T2LNum(cnum)};
                    tdata=[tdata; pos];
                    tdata=removetabel(tdata);
                    set(htable, 'Data', tdata); 
                    set(data.chls(cnum), 'LineWidth', 1.5, 'Color', 'm');
                    pause(0.5); 
                    set(data.chls(cnum), 'LineWidth', 0.5, 'Color', 'k');        
                end
            end
        case 'Mirror'
            htable=argin{2};
            hedit=argin{4};
            
            empt=str2num(get(hedit, 'String'));  
            cnum=str2num(get(callbackdata.Source, 'Tag'));
            data=get(handles.axes1,'UserData');            
            tdata=htable.Data;
            if strcmp(data.crvs(cnum).form, 'L-NURL')
                if isempty(empt)
                    lnum=data.T2LNum(cnum);
                    set(hedit, 'String', num2str(lnum));
                    set(data.chls(cnum), 'LineWidth', 1.5, 'Color', 'm');
                    pause(0.5); 
                    set(data.chls(cnum), 'LineWidth', 0.5, 'Color', 'k'); 
                else
                    pos={'Curve', data.T2LNum(cnum)};
                    tdata=[tdata; pos];
                    tdata=removetabel(tdata);
                    set(htable, 'Data', tdata); 
                    set(data.chls(cnum), 'LineWidth', 1.5, 'Color', 'm');
                    pause(0.5); 
                    set(data.chls(cnum), 'LineWidth', 0.5, 'Color', 'k');        
                end
            end
        case {'Split', 'Glue'}
            htable=argin{2};
            
            cnum=str2num(get(callbackdata.Source, 'Tag'));
            data=get(handles.axes1,'UserData');            
            tdata=htable.Data;
            if strcmp(data.crvs(cnum).form, 'L-NURL')
                pos={'Curve', data.T2LNum(cnum)};
                tdata=[tdata; pos];
                tdata=removetabel(tdata);
                set(htable, 'Data', tdata); 
                set(data.chls(cnum), 'LineWidth', 1.5, 'Color', 'm');
                pause(0.5); 
                set(data.chls(cnum), 'LineWidth', 0.5, 'Color', 'k');        
            end
        case {'Translate', 'Rotate', 'Double', 'Save', 'Scale'}
            htable=argin{2};
            
            cnum=str2num(get(callbackdata.Source, 'Tag'));
            data=get(handles.axes1,'UserData');            
            tdata=htable.Data;
            if strcmp(data.crvs(cnum).form, 'L-NURL')
                pos={'Curve', data.T2LNum(cnum)};
                tdata=[tdata; pos];
                tdata=removetabel(tdata);
                set(htable, 'Data', tdata); 
                set(data.chls(cnum), 'LineWidth', 1.5, 'Color', 'm');
                pause(0.5); 
                set(data.chls(cnum), 'LineWidth', 0.5, 'Color', 'k');        
            end
        case 'CutByPoint'
            htable=argin{2};
            hedit=argin{4};
            
            Type=get(callbackdata.Source, 'Type');
            cnum=str2num(get(callbackdata.Source, 'Tag'));
            data=get(handles.axes1,'UserData');
            tdata=htable.Data;
            if strcmp(data.crvs(cnum).form, 'POINTS') || strcmp(Type, 'text')
                pnum=data.T2PNum(cnum); 
                set(hedit, 'String', num2str(pnum));
                set(data.chls(cnum), 'LineWidth', 1.5, 'Color', 'm');
                pause(0.5); 
                set(data.chls(cnum), 'LineWidth', 0.5, 'Color', 'k'); 
            elseif strcmp(data.crvs(cnum).form, 'L-NURL')
                pos={'Curve', data.T2LNum(cnum)};
                n=numel(tdata)/2; 
                if n<1
                    tdata=[tdata; pos];
                    set(htable, 'Data', tdata); 
                end
                set(data.chls(cnum), 'LineWidth', 1.5, 'Color', 'm');
                pause(0.5); 
                set(data.chls(cnum), 'LineWidth', 0.5, 'Color', 'k');        
            end        
        case {'Extract', 'Center'}
            Type=get(callbackdata.Source, 'Type');
            cnum=str2num(get(callbackdata.Source, 'Tag'));
            data=get(handles.axes1,'UserData');
            htable=argin{2};
            tdata=htable.Data;
            if strcmp(Type, 'line')
                pos={'Curve', data.T2LNum(cnum)};
                tdata=[tdata; pos];
                tdata=removetabel(tdata);
                set(htable, 'Data', tdata); 
                set(data.chls(cnum), 'LineWidth', 1.5, 'Color', 'm');
                pause(0.5); 
                set(data.chls(cnum), 'LineWidth', 0.5, 'Color', 'k');                
            end
        case {'LinePN', 'LinePT'}
            hedit=argin{2};
            
            Type=get(callbackdata.Source, 'Type');
            cnum=str2num(get(callbackdata.Source, 'Tag'));
            data=get(handles.axes1,'UserData');    
            if strcmp(Type, 'line')
                lnum=data.T2LNum(cnum);
                set(hedit, 'String', num2str(lnum));
                set(data.chls(cnum), 'LineWidth', 1.5, 'Color', 'm');
                pause(0.5); 
                set(data.chls(cnum), 'LineWidth', 0.5, 'Color', 'k'); 
            end
    end

    
% Remove data from a table of two columns
function rdata=removetabel(tdata)
n=numel(tdata)/2;
Num=zeros(n,1);
for i=1:n
    Num(i)=tdata{n+i};
end
[~, u]=RemDuplicate(Num);
u=find(u);
uu=[u,n+u];
n=length(uu)/2;
rdata=cell(n,2);
for i=1:n
    rdata{i}=tdata{uu(i)};
    rdata{n+i}=tdata{uu(n+i)};
end

    
% Edit curves
function Edit(hObject, callbackdata, handles, argin) 
    switch argin{1}
        case 'Deleted'
            htable=argin{2};
            tdata=htable.Data;
            n=numel(tdata)/2;
            data=get(handles.axes1,'UserData'); 
            if n>0
                Type=cell(1, n);
                Num=cell(1, n);
                for i=1:n
                    Type{i}=tdata{i};
                    Num{i}=tdata{n+i};
                end                
                r=1; s=1; t=1;
                T=[]; dhL=[]; dhP=[];
                for i=1:n
                    if strcmp(Type{i}, 'Curve')
                        T(r)=data.LNumber(Num{i});
                        dhL(s)=Num{i};
                        r=r+1; s=s+1;
                    elseif strcmp(Type{i}, 'Point')
                        T(r)=data.PNumber(Num{i});
                        dhP(t)=Num{i};
                        r=r+1; t=t+1;
                    end
                end
                delete(data.chls(T)); 
                data.dcrvs=[data.dcrvs, data.crvs(T)]; 
                data.disline=[data.disline, data.isline(T)]; 
                data.crvs(T)=[]; 
                data.chls(T)=[]; 
                data.isline(T)=[]; 
                set(handles.axes1, 'UserData', data); 
            end
            
            % Empty the list of deleted curves
            Reset(hObject, callbackdata, handles, argin);
            
            % Refresh line and point numbers
            RefreshPlot_Callback(hObject, callbackdata, handles);
            
            % Assign function handle to point numbers and curve handles
            data=get(handles.axes1,'UserData'); 
            set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, argin});
            set(handles.axes1,'UserData', data);       
        case 'Save'
            htable=argin{2};
            tdata=htable.Data;
            hedit=argin{4};
            n=numel(tdata)/2;
            data=get(handles.axes1,'UserData'); 
            if n>0
                Type=cell(1, n);
                Num=cell(1, n);
                for i=1:n
                    Type{i}=tdata{i};
                    Num{i}=tdata{n+i};
                end                
                r=1; s=1; t=1;
                T=[]; dhL=[]; dhP=[];
                for i=1:n
                    if strcmp(Type{i}, 'Curve')
                        T(r)=data.LNumber(Num{i});
                        dhL(s)=Num{i};
                        r=r+1; s=s+1;
                    elseif strcmp(Type{i}, 'Point')
                        T(r)=data.PNumber(Num{i});
                        dhP(t)=Num{i};
                        r=r+1; t=t+1;
                    end
                end
                sdata.crvs=data.crvs(T);
                sdata.isline=data.isline(T);
                filename=get(hedit, 'String');
                save(filename, 'sdata');
            end
            
            % Empty the list of deleted curves
            Reset(hObject, callbackdata, handles, argin);            
        case 'Fillet'
            htable=argin{2}; 
            hedit=argin{4};
            
            tdata=htable.Data; 
            n=numel(tdata)/2; 
            data=get(handles.axes1,'UserData'); 
            if n==2
                % Delete the two lines from data
                Num=cell(1, n);
                for i=1:n
                    Num{i}=tdata{n+i};
                end                
                r=1;  T=[]; dhL=[]; 
                for i=1:n
                    T(r)=data.LNumber(Num{i});
                    dhL(r)=Num{i};
                    r=r+1;
                end
                lines=data.crvs(T);                
                
                % Fillet the two lines to get three new lines2points
                R=str2num(get(hedit, 'String')); 
                [crvs, center]=nrlfillet(lines(1), lines(2), R);
                cpnt=pntmake(center');
                phls=pntplot(center, 0.5, '.');
                isline=[1, 1, 1, 0];
                chls = nrlcrvplot(crvs,50);
                crvs=[crvs, cpnt];
                chls=[chls, phls];
                
                % Set the axes limite
                XLim=get(handles.axes1,'XLim');
                YLim=get(handles.axes1,'YLim');
                set(handles.editXLim,'String', num2str(XLim));
                set(handles.editYLim,'String', num2str(YLim));
                
                % Select curves for further operation
                set(chls, 'ButtonDownFcn', {@Selected, handles, argin}); 
                
                % Set Tag to line handls
                m=length(data.crvs);
                for i=1:3
                    set(chls(i), 'Tag', num2str(m+i));
                end
                
                % Delete the used lines
                delete(data.chls(T)); 
                data.dcrvs=[data.dcrvs, data.crvs(T)]; 
                data.disline=[data.disline, data.isline(T)]; 
                data.crvs(T)=[]; 
                data.chls(T)=[]; 
                data.isline(T)=[]; 
                
                % Save data to 'UserData' of handles.axes1
                data.crvs=[data.crvs, crvs];  
                data.chls=[data.chls, chls]; 
                data.isline=[data.isline, isline]; 
                data.n=data.n+1;
                set(handles.axes1, 'UserData', data);  
                
                % Refresh line and point numbers
                RefreshPlot_Callback(hObject, callbackdata, handles);
            end
            
            % Empty the list of deleted curves
            Reset(hObject, callbackdata, handles, argin);     
        case 'Measure'
            htable=argin{2}; 
            hedit=argin{4};
            
            tdata=htable.Data; 
            n=numel(tdata)/2; 
            data=get(handles.axes1,'UserData'); 
            if n==1
                % Get the curve from data
                Num=tdata{2};  
                T=data.LNumber(Num);
                crv=data.crvs(T);     
                dist= nrlmeasure (crv);
                
                % Fillet the two lines2points to get three new lines2points
                set(hedit, 'String', num2str(dist)); 
            end            
        case 'FilletTA'
            htable=argin{2}; 
            hedit=argin{4};
            
            tdata=htable.Data; 
            n=numel(tdata)/2; 
            data=get(handles.axes1,'UserData'); 
            S=str2num(get(hedit, 'String')); 
            if n==2 && ~isempty(S)
                % Delete the two lines2points from data
                Num=cell(1, n);
                for i=1:n
                    Num{i}=tdata{n+i};
                end                
                r=1;  T=[]; 
                for i=1:n
                    T(r)=data.LNumber(Num{i});
                    r=r+1;
                end
                lines=data.crvs(T);                
                
                % Fillet the two lines2points to get three new lines2points
                n=length(S);
                if n==1
                    [crvs, center]=nrlfillet(lines(1), lines(2), S);
                    cpnt=pntmake(center');
                    phls=pntplot(center, 0.5, '.');
                    isline=[1, 1, 1, 0];
                    chls = nrlcrvplot(crvs,50);
                    crvs=[crvs, cpnt];
                    chls=[chls, phls];
                elseif S(1) && ~S(2)
                    [crvs, center]=tangentarcs(lines(1), lines(2));
                    cpnt=pntmake(center);
                    phls=pntplot(center, 0.5, '.');
                    t=length(crvs);
                    isline=[ones(1, t), 0];
                    chls = nrlcrvplot(crvs,50);
                    crvs=[crvs, cpnt];
                    chls=[chls, phls];
                elseif ~S(1) && S(2)
                    [~, ~, crvs, center]=tangentarcs(lines(1), lines(2));
                    cpnt=pntmake(center);
                    phls=pntplot(center, 0.5, '.');
                    t=length(crvs);
                    isline=[ones(1, t), 0];
                    chls = nrlcrvplot(crvs,50);
                    crvs=[crvs, cpnt];
                    chls=[chls, phls];
                elseif S(1) && S(2)
                    [~, center1, ~, center2, crvs]=tangentarcs(lines(1), lines(2));
                    cpnt1=pntmake(center1);
                    phls1=pntplot(center1, 0.5, '.');
                    cpnt2=pntmake(center2);
                    phls2=pntplot(center2, 0.5, '.');
                    t=length(crvs);
                    isline=[ones(1, t), 0, 0];
                    chls = nrlcrvplot(crvs,50);
                    crvs=[crvs, cpnt1, cpnt2];
                    chls=[chls, phls1, phls2];
                end
                
                % Set the axes limite
                XLim=get(handles.axes1,'XLim');
                YLim=get(handles.axes1,'YLim');
                set(handles.editXLim,'String', num2str(XLim));
                set(handles.editYLim,'String', num2str(YLim));
                
                % Select curves for further operation
                set(chls, 'ButtonDownFcn', {@Selected, handles, argin}); 
                
                % Set Tag to line handls
                m=length(data.crvs);
                s=length(crvs);
                for i=1:s
                    set(chls(i), 'Tag', num2str(m+i));
                end
                
                % Delete the used lines2points
                delete(data.chls(T)); 
                data.dcrvs=[data.dcrvs, data.crvs(T)]; 
                data.disline=[data.disline, data.isline(T)]; 
                data.crvs(T)=[]; 
                data.chls(T)=[]; 
                data.isline(T)=[]; 
                
                % Save data to 'UserData' of handles.axes1
                data.crvs=[data.crvs, crvs];  
                data.chls=[data.chls, chls]; 
                data.isline=[data.isline, isline]; 
                data.n=data.n+1;
                set(handles.axes1, 'UserData', data);  
                
                % Refresh line and point numbers
                RefreshPlot_Callback(hObject, callbackdata, handles);
            end
            
            % Empty the list of deleted curves
            Reset(hObject, callbackdata, handles, argin);     
        case 'Cut'
            htable=argin{2}; 
            hedit=argin{4};     
            
            % Get the curve to cut
            data=get(handles.axes1,'UserData'); 
            lnum=str2num(get(hedit, 'String')); 
            t=data.LNumber(lnum);
            crv=data.crvs(t);
            
            % Get the curves to be cutted
            tdata=htable.Data; 
            n=numel(tdata)/2; 
            if n>=1
                % Delete the curves from data
                Num=cell(1, n);
                for i=1:n
                    Num{i}=tdata{n+i};
                end
                r=1;  T=[]; 
                for i=1:n
                    T(r)=data.LNumber(Num{i});
                    r=r+1;
                end
                crvs=data.crvs(T);
                T=[T, t];
                                
                % Cut the curves to get new curves
                [crvn, crvsn]=nrlcuts(crv, crvs);
                crvs=[crvn, crvsn];
                n=numel(crvs);
                isline=ones(1,n);
                chls = nrlcrvplot(crvs,50);
                
                % Set the axes limite
                XLim=get(handles.axes1,'XLim');
                YLim=get(handles.axes1,'YLim');
                set(handles.editXLim,'String', num2str(XLim));
                set(handles.editYLim,'String', num2str(YLim));
                
                % Select curves for further operation
                set(chls, 'ButtonDownFcn', {@Selected, handles, argin}); 
                
                % Set Tag to line handls
                m=length(data.crvs);
                for i=1:n
                    set(chls(i), 'Tag', num2str(m+i));
                end
                
                % Delete the used curves
                delete(data.chls(T)); 
                data.dcrvs=[data.dcrvs, data.crvs(T)]; 
                data.disline=[data.disline, data.isline(T)]; 
                data.crvs(T)=[]; 
                data.chls(T)=[]; 
                data.isline(T)=[]; 
                
                % Save data to 'UserData' of handles.axes1
                data.crvs=[data.crvs, crvs];  
                data.chls=[data.chls, chls]; 
                data.isline=[data.isline, isline]; 
                data.n=data.n+1;
                set(handles.axes1, 'UserData', data);  
                
                % Refresh line and point numbers
                RefreshPlot_Callback(hObject, callbackdata, handles);
            end
            
            % Empty the list of deleted curves
            Reset(hObject, callbackdata, handles, argin);     
        case 'Mirror'
            htable=argin{2}; 
            hedit=argin{4};     
            
            % Get the reference curve
            data=get(handles.axes1,'UserData'); 
            lnum=str2num(get(hedit, 'String')); 
            t=data.LNumber(lnum);
            crv=data.crvs(t);
            
            % Get the curves to be mirrored
            tdata=htable.Data; 
            n=numel(tdata)/2; 
            if n>=1
                % Get the curves from data
                Num=cell(1, n);
                for i=1:n
                    Num{i}=tdata{n+i};
                end
                r=1;  T=[]; 
                for i=1:n
                    T(r)=data.LNumber(Num{i});
                    r=r+1;
                end
                crvs=data.crvs(T);
                
                % Mirror the curves to get new curves
                crvs=plcrvmirror(crvs, crv);
                n=numel(crvs);
                isline=ones(1,n);
                chls = nrlcrvplot(crvs,50);
                axis equal;
                
                % Set the axes limite
                XLim=get(handles.axes1,'XLim');
                YLim=get(handles.axes1,'YLim');
                set(handles.editXLim,'String', num2str(XLim));
                set(handles.editYLim,'String', num2str(YLim));
                
                % Select curves for further operation
                set(chls, 'ButtonDownFcn', {@Selected, handles, argin}); 
                
                % Set Tag to line handls
                m=length(data.crvs);
                for i=1:n
                    set(chls(i), 'Tag', num2str(m+i));
                end
                
                % Save data to 'UserData' of handles.axes1
                data.crvs=[data.crvs, crvs];  
                data.chls=[data.chls, chls]; 
                data.isline=[data.isline, isline]; 
                data.n=data.n+1;
                set(handles.axes1, 'UserData', data);  
                
                % Refresh line and point numbers
                RefreshPlot_Callback(hObject, callbackdata, handles);
            end
            
            % Empty the list of deleted curves
            Reset(hObject, callbackdata, handles, argin);     
        case 'Split'
            htable=argin{2}; 
            hedit=argin{4};     
            
            % Get the knots to split
            data=get(handles.axes1,'UserData'); 
            t=str2num(get(hedit, 'String'));             
            if length(t)==1
                if t>=2
                    t=linspace(0, 1, fix(t+1));
                end
            end
            
            % Get the curves to be splitted
            tdata=htable.Data; 
            n=numel(tdata)/2; 
            if n>=1
                % Delete the curves from data
                Num=cell(1, n);
                for i=1:n
                    Num{i}=tdata{n+i};
                end
                r=1;  T=[]; 
                for i=1:n
                    T(r)=data.LNumber(Num{i});
                    r=r+1;
                end
                crvns=data.crvs(T);
                                
                % Split the curves to get new curves
                crvs=[];
                for i=1:n
                    crvis=nrlsplits(crvns(i), t);
                    crvs=[crvs, crvis];
                end        
                n=numel(crvs);
                isline=ones(1,n);
                chls = nrlcrvplot(crvs,50);
                
                % Set the axes limite
                XLim=get(handles.axes1,'XLim');
                YLim=get(handles.axes1,'YLim');
                set(handles.editXLim,'String', num2str(XLim));
                set(handles.editYLim,'String', num2str(YLim));
                
                % Select curves for further operation
                set(chls, 'ButtonDownFcn', {@Selected, handles, argin}); 
                
                % Set Tag to line handls
                m=length(data.crvs);
                for i=1:n
                    set(chls(i), 'Tag', num2str(m+i));
                end
                
                % Delete the used curves
                delete(data.chls(T)); 
                data.dcrvs=[data.dcrvs, data.crvs(T)]; 
                data.disline=[data.disline, data.isline(T)]; 
                data.crvs(T)=[]; 
                data.chls(T)=[]; 
                data.isline(T)=[]; 
                
                % Save data to 'UserData' of handles.axes1
                data.crvs=[data.crvs, crvs];  
                data.chls=[data.chls, chls]; 
                data.isline=[data.isline, isline]; 
                data.n=data.n+1;
                set(handles.axes1, 'UserData', data);  
                
                % Refresh line and point numbers
                RefreshPlot_Callback(hObject, callbackdata, handles);
            end
            
            % Empty the list of deleted curves
            Reset(hObject, callbackdata, handles, argin);    
        case 'Translate'
            htable=argin{2}; 
            hedit=argin{4};     
            hedit2=argin{5};  
            
            % Get the plane vector to trans late
            data=get(handles.axes1,'UserData'); 
            vec=str2num(get(hedit, 'String'));  
            isrt=str2num(get(hedit2, 'String'));  
            
            % Get the curves to be splitted
            tdata=htable.Data; 
            n=numel(tdata)/2; 
            if n>=1
                % Delete the curves from data
                Num=cell(1, n);
                for i=1:n
                    Num{i}=tdata{n+i};
                end
                r=1;  T=[]; 
                for i=1:n
                    T(r)=data.LNumber(Num{i});
                    r=r+1;
                end
                crvns=data.crvs(T);
                                
                % Split the curves to get new curves
                crvs=[];
                for i=1:n
                    crvis=plcrvtranslat(crvns(i), vec);
                    crvs=[crvs, crvis];
                end
                n=numel(crvs);
                isline=ones(1,n);
                chls = nrlcrvplot(crvs,50);
                
                % Set the axes limite
                XLim=get(handles.axes1,'XLim');
                YLim=get(handles.axes1,'YLim');
                set(handles.editXLim,'String', num2str(XLim));
                set(handles.editYLim,'String', num2str(YLim));
                
                % Select curves for further operation
                set(chls, 'ButtonDownFcn', {@Selected, handles, argin}); 
                
                % Set Tag to line handls
                m=length(data.crvs);
                for i=1:n
                    set(chls(i), 'Tag', num2str(m+i));
                end
                
                % Delete the used curves
                if isrt==0
                    delete(data.chls(T)); 
                    data.dcrvs=[data.dcrvs, data.crvs(T)]; 
                    data.disline=[data.disline, data.isline(T)]; 
                    data.crvs(T)=[]; 
                    data.chls(T)=[]; 
                    data.isline(T)=[]; 
                end
                
                % Save data to 'UserData' of handles.axes1
                data.crvs=[data.crvs, crvs];  
                data.chls=[data.chls, chls]; 
                data.isline=[data.isline, isline]; 
                data.n=data.n+1;
                set(handles.axes1, 'UserData', data);  
                
                % Refresh line and point numbers
                RefreshPlot_Callback(hObject, callbackdata, handles);
            end
            
            % Empty the list of deleted curves
            Reset(hObject, callbackdata, handles, argin);  
        case 'Scale'
            htable=argin{2}; 
            hedit=argin{4};     
            hedit2=argin{5};  
            
            % Get the plane vector to scale
            data=get(handles.axes1,'UserData'); 
            scale=str2num(get(hedit, 'String'));  
            isrt=str2num(get(hedit2, 'String'));  
            vec=[scale, 1];
            sm= vecscale(vec);
            
            % Get the curves to be splitted
            tdata=htable.Data; 
            n=numel(tdata)/2; 
            if n>=1
                % Delete the curves from data
                Num=cell(1, n);
                for i=1:n
                    Num{i}=tdata{n+i};
                end
                r=1;  T=[]; 
                for i=1:n
                    T(r)=data.LNumber(Num{i});
                    r=r+1;
                end
                crvns=data.crvs(T);
                                
                % Split the curves to get new curves
                crvs=[];
                for i=1:n
                    crvis=nrltform(crvns(i), sm);
                    crvs=[crvs, crvis];
                end
                n=numel(crvs);
                isline=ones(1,n);
                chls = nrlcrvplot(crvs,50);
                
                % Set the axes limite
                XLim=get(handles.axes1,'XLim');
                YLim=get(handles.axes1,'YLim');
                set(handles.editXLim,'String', num2str(XLim));
                set(handles.editYLim,'String', num2str(YLim));
                
                % Select curves for further operation
                set(chls, 'ButtonDownFcn', {@Selected, handles, argin}); 
                
                % Set Tag to line handls
                m=length(data.crvs);
                for i=1:n
                    set(chls(i), 'Tag', num2str(m+i));
                end
                
                % Delete the used curves
                if isrt==0
                    delete(data.chls(T)); 
                    data.dcrvs=[data.dcrvs, data.crvs(T)]; 
                    data.disline=[data.disline, data.isline(T)]; 
                    data.crvs(T)=[]; 
                    data.chls(T)=[]; 
                    data.isline(T)=[]; 
                end
                
                % Save data to 'UserData' of handles.axes1
                data.crvs=[data.crvs, crvs];  
                data.chls=[data.chls, chls]; 
                data.isline=[data.isline, isline]; 
                data.n=data.n+1;
                set(handles.axes1, 'UserData', data);  
                
                % Refresh line and point numbers
                RefreshPlot_Callback(hObject, callbackdata, handles);
            end
            
            % Empty the list of deleted curves
            Reset(hObject, callbackdata, handles, argin);  
        case 'Double'
            htable=argin{2}; 
            hedit=argin{4};     
            hedit2=argin{5};  
            
            % Get the plane vector to trans late
            data=get(handles.axes1,'UserData'); 
            np=str2num(get(hedit, 'String'));  
            dt=str2num(get(hedit2, 'String'));  
            
            % Get the curves to be splitted
            tdata=htable.Data; 
            n=numel(tdata)/2; 
            if n>=1
                % Delete the curves from data
                Num=cell(1, n);
                for i=1:n
                    Num{i}=tdata{n+i};
                end
                r=1;  T=[]; 
                for i=1:n
                    T(r)=data.LNumber(Num{i});
                    r=r+1;
                end
                crvns=data.crvs(T);
                                
                % Double the curves to get new curves
                crvs=[];
                for i=1:n
                    [crv1, crv2]=bertrand(crvns(i), dt, np);
                    crvs=[crvs, crv1, crv2];
                end
                n=numel(crvs);
                isline=ones(1,n);
                chls = nrlcrvplot(crvs,50);
                
                % Set the axes limite
                XLim=get(handles.axes1,'XLim');
                YLim=get(handles.axes1,'YLim');
                set(handles.editXLim,'String', num2str(XLim));
                set(handles.editYLim,'String', num2str(YLim));
                
                % Select curves for further operation
                set(chls, 'ButtonDownFcn', {@Selected, handles, argin}); 
                
                % Set Tag to line handls
                m=length(data.crvs);
                for i=1:n
                    set(chls(i), 'Tag', num2str(m+i));
                end
                                
                % Save data to 'UserData' of handles.axes1
                data.crvs=[data.crvs, crvs];  
                data.chls=[data.chls, chls]; 
                data.isline=[data.isline, isline]; 
                data.n=data.n+1;
                set(handles.axes1, 'UserData', data);  
                
                % Refresh line and point numbers
                RefreshPlot_Callback(hObject, callbackdata, handles);
            end
            
            % Empty the list of deleted curves
            Reset(hObject, callbackdata, handles, argin);    
        case 'Rotate'
            htable=argin{2}; 
            hedit=argin{4};     
            hedit2=argin{5};  
            hedit3=argin{6};  
            
            % Get the plane vector to trans late
            data=get(handles.axes1,'UserData'); 
            angle=str2num(get(hedit, 'String'));  
            pnt=str2num(get(hedit2, 'String'));  
            isrt=str2num(get(hedit3, 'String'));  
            
            % Get the curves to be splitted
            tdata=htable.Data; 
            n=numel(tdata)/2; 
            if n>=1
                % Delete the curves from data
                Num=cell(1, n);
                for i=1:n
                    Num{i}=tdata{n+i};
                end
                r=1;  T=[]; 
                for i=1:n
                    T(r)=data.LNumber(Num{i});
                    r=r+1;
                end
                crvns=data.crvs(T);
                                
                % Split the curves to get new curves
                pnt=GetCoords(handles, pnt);
                crvs=[];
                for i=1:n
                    crvis=crvrotate(crvns(i), angle, pnt);
                    crvs=[crvs, crvis];
                end
                n=numel(crvs);
                isline=ones(1,n);
                chls = nrlcrvplot(crvs,50);
                
                % Set the axes limite
                XLim=get(handles.axes1,'XLim');
                YLim=get(handles.axes1,'YLim');
                set(handles.editXLim,'String', num2str(XLim));
                set(handles.editYLim,'String', num2str(YLim));
                
                % Select curves for further operation
                set(chls, 'ButtonDownFcn', {@Selected, handles, argin}); 
                
                % Set Tag to line handls
                m=length(data.crvs);
                for i=1:n
                    set(chls(i), 'Tag', num2str(m+i));
                end
                
                % Delete the used curves
                if isrt==0
                    delete(data.chls(T)); 
                    data.dcrvs=[data.dcrvs, data.crvs(T)]; 
                    data.disline=[data.disline, data.isline(T)]; 
                    data.crvs(T)=[]; 
                    data.chls(T)=[]; 
                    data.isline(T)=[]; 
                end
                
                % Save data to 'UserData' of handles.axes1
                data.crvs=[data.crvs, crvs];  
                data.chls=[data.chls, chls]; 
                data.isline=[data.isline, isline]; 
                data.n=data.n+1;
                set(handles.axes1, 'UserData', data);  
                
                % Refresh line and point numbers
                RefreshPlot_Callback(hObject, callbackdata, handles);
            end
            
            % Empty the list of deleted curves
            Reset(hObject, callbackdata, handles, argin);    
        case 'Glue'
            htable=argin{2}; 
            
            % Get the curves to be glued
            data=get(handles.axes1,'UserData'); 
            tdata=htable.Data; 
            n=numel(tdata)/2; 
            if n>=1
                % Delete the curves from data
                Num=cell(1, n);
                for i=1:n
                    Num{i}=tdata{n+i};
                end
                r=1;  T=[]; 
                for i=1:n
                    T(r)=data.LNumber(Num{i});
                    r=r+1;
                end
                crvns=data.crvs(T);
                                
                % Glue the curves to get a new curve
                crvs=nrlglues(crvns);
                n=numel(crvs);
                isline=ones(1,n);
                chls = nrlcrvplot(crvs,50);
                
                % Set the axes limite
                XLim=get(handles.axes1,'XLim');
                YLim=get(handles.axes1,'YLim');
                set(handles.editXLim,'String', num2str(XLim));
                set(handles.editYLim,'String', num2str(YLim));
                
                % Select curves for further operation
                set(chls, 'ButtonDownFcn', {@Selected, handles, argin}); 
                
                % Set Tag to line handls
                m=length(data.crvs);
                for i=1:n
                    set(chls(i), 'Tag', num2str(m+i));
                end
                
                % Delete the used curves
                delete(data.chls(T)); 
                data.dcrvs=[data.dcrvs, data.crvs(T)]; 
                data.disline=[data.disline, data.isline(T)]; 
                data.crvs(T)=[]; 
                data.chls(T)=[]; 
                data.isline(T)=[]; 
                
                % Save data to 'UserData' of handles.axes1
                data.crvs=[data.crvs, crvs];  
                data.chls=[data.chls, chls]; 
                data.isline=[data.isline, isline]; 
                data.n=data.n+1;
                set(handles.axes1, 'UserData', data);  
                
                % Refresh line and point numbers
                RefreshPlot_Callback(hObject, callbackdata, handles);
            end
            
            % Empty the list of deleted curves
            Reset(hObject, callbackdata, handles, argin);    
        case 'CutByPoint'
            htable=argin{2}; 
            hedit=argin{4};     
            
            % Get the curve to cut
            data=get(handles.axes1,'UserData'); 
            pnum=str2num(get(hedit, 'String')); 
            t=data.PNumber(pnum);
            pnt=data.crvs(t).coefs;
            
            % Get the curves to be cutted
            tdata=htable.Data;
            n=numel(tdata)/2;
            if n>=1
                % Delete the curve from data
                Num=tdata{n+1};
                T=data.LNumber(Num);
                crv=data.crvs(T);
                
                % Get the nearest point from the curve
                [~, um, ~]=crvnearpnt(crv, pnt');
                
                % Split the curve by the nearest point of the curve
                [crv1, crv2]=nrlsplit(crv, um);
                crvs=[crv1, crv2];
                n=numel(crvs);
                isline=ones(1,n);
                chls = nrlcrvplot(crvs,50);
                
                % Set the axes limite
                XLim=get(handles.axes1,'XLim');
                YLim=get(handles.axes1,'YLim');
                set(handles.editXLim,'String', num2str(XLim));
                set(handles.editYLim,'String', num2str(YLim));
                
                % Select curves for further operation
                set(chls, 'ButtonDownFcn', {@Selected, handles, argin}); 
                
                % Set Tag to line handls
                m=length(data.crvs);
                for i=1:n
                    set(chls(i), 'Tag', num2str(m+i));
                end
                
                % Delete the used curve
                delete(data.chls(T)); 
                data.dcrvs=[data.dcrvs, data.crvs(T)]; 
                data.disline=[data.disline, data.isline(T)]; 
                data.crvs(T)=[]; 
                data.chls(T)=[]; 
                data.isline(T)=[]; 
                
                % Save data to 'UserData' of handles.axes1
                data.crvs=[data.crvs, crvs];  
                data.chls=[data.chls, chls]; 
                data.isline=[data.isline, isline]; 
                data.n=data.n+1;
                set(handles.axes1, 'UserData', data);  
                
                % Refresh line and point numbers
                RefreshPlot_Callback(hObject, callbackdata, handles);
            end
            
            % Empty the list of deleted curves
            Reset(hObject, callbackdata, handles, argin);       
        case 'Extract'
            htable=argin{2}; 
            
            tdata=htable.Data; 
            n=numel(tdata)/2; 
            data=get(handles.axes1,'UserData'); 
            if n>0
                % Find the curves from data
                Num=cell(1, n);
                for i=1:n
                    Num{i}=tdata{n+i};
                end                
                r=1;  T=[]; 
                for i=1:n
                    T(r)=data.LNumber(Num{i});
                    r=r+1;
                end
                crvn=data.crvs(T);                
                
                % Fillet the two lines2points to get three new lines2points
                crvs=[]; chls=[];
                for i=1:n
                    pnts=nrleval(crvn(i), [0, 0.5, 1]);
                    cpnt1=pntmake(pnts(1:2,1)');
                    phl1=pntplot(pnts(1:2,1), 0.5, '.');
                    cpnt2=pntmake(pnts(1:2,2)');
                    phl2=pntplot(pnts(1:2,2), 0.5, '.');
                    cpnt3=pntmake(pnts(1:2,3)');
                    phl3=pntplot(pnts(1:2,3), 0.5, '.');
                    crvs=[crvs, cpnt1, cpnt2, cpnt3];
                    chls=[chls, phl1, phl2, phl3];
                end
                isline=zeros(1,3*n);
                
                % Set the axes limite
                XLim=get(handles.axes1,'XLim');
                YLim=get(handles.axes1,'YLim');
                set(handles.editXLim,'String', num2str(XLim));
                set(handles.editYLim,'String', num2str(YLim));
                
                % Select curves for further operation
                set(chls, 'ButtonDownFcn', {@Selected, handles, argin}); 
                
                % Set Tag to line handls
                m=length(data.crvs);
                for i=1:3*n
                    set(chls(i), 'Tag', num2str(m+i));
                end
                
                % Save data to 'UserData' of handles.axes1
                data.crvs=[data.crvs, crvs];  
                data.chls=[data.chls, chls]; 
                data.isline=[data.isline, isline]; 
                data.n=data.n+1;
                set(handles.axes1, 'UserData', data);  
                
                % Refresh line and point numbers
                RefreshPlot_Callback(hObject, callbackdata, handles);
                
                % Remove dumplicated points
                ReDpPoints(hObject, callbackdata, handles, argin);
            end
            
            % Empty the list of deleted curves
            Reset(hObject, callbackdata, handles, argin);    
        case 'Center'
            htable=argin{2}; 
            
            tdata=htable.Data; 
            n=numel(tdata)/2; 
            data=get(handles.axes1,'UserData'); 
            if n>0
                % Find the curves from data
                Num=cell(1, n);
                for i=1:n
                    Num{i}=tdata{n+i};
                end                
                r=1;  T=[]; 
                for i=1:n
                    T(r)=data.LNumber(Num{i});
                    r=r+1;
                end
                crvn=data.crvs(T);                
                
                % Fillet the two lines2points to get three new lines2points
                crvs=[]; chls=[];
                for i=1:n
                    center = nrlcircenter(crvn(i));
                    cpnt=pntmake(center);
                    phl=pntplot(center, 0.5, '.');
                    crvs=[crvs, cpnt];
                    chls=[chls, phl];
                end
                isline=zeros(1,n);
                
                % Set the axes limite
                XLim=get(handles.axes1,'XLim');
                YLim=get(handles.axes1,'YLim');
                set(handles.editXLim,'String', num2str(XLim));
                set(handles.editYLim,'String', num2str(YLim));
                
                % Select curves for further operation
                set(chls, 'ButtonDownFcn', {@Selected, handles, argin}); 
                
                % Set Tag to line handls
                m=length(data.crvs);
                for i=1:n
                    set(chls(i), 'Tag', num2str(m+i));
                end
                
                % Save data to 'UserData' of handles.axes1
                data.crvs=[data.crvs, crvs];  
                data.chls=[data.chls, chls]; 
                data.isline=[data.isline, isline]; 
                data.n=data.n+1;
                set(handles.axes1, 'UserData', data);  
                
                % Refresh line and point numbers
                RefreshPlot_Callback(hObject, callbackdata, handles);
            end
            
            % Empty the list of deleted curves
            Reset(hObject, callbackdata, handles, argin);    
        case 'LinePN'
            hedit2=argin{2}; 
            hedit1=argin{4};     
            hedit3=argin{5};     
            
            % Get the point and the length of the line
            pnt1=str2num(get(hedit1, 'String'));
            pnt1=GetCoords(handles, pnt1);
            Len=str2num(get(hedit3, 'String'));
            
            % Get the normal line
            data=get(handles.axes1,'UserData'); 
            lnum=str2num(get(hedit2, 'String')); 
            t=data.LNumber(lnum);
            crv=data.crvs(t);
            
            % Get the normal at the point
            pnt1=[pnt1, 0]';
            dr=neartangent(crv, pnt1);
            pnt2=pnt1+Len*dr;
            crvs=nrlline(pnt1', pnt2');
            chls = nrlcrvplot(crvs, 50);  

            % Set the axes limite
            XLim=get(handles.axes1,'XLim');
            YLim=get(handles.axes1,'YLim');
            set(handles.editXLim,'String', num2str(XLim));
            set(handles.editYLim,'String', num2str(YLim));

            % Select curves for further operation
            set(chls, 'ButtonDownFcn', {@Selected, handles, {'Selected'}});

            % Judge whether the entity is a line or a point
            n=length(chls);
            isline=ones(1,n);

            % Set Tag to line handls
            m=length(data.crvs);
            for i=1:n
                set(chls(i), 'Tag', num2str(m+i));
            end

            % Save data to 'UserData' of handles.axes1
            data.crvs=[data.crvs, crvs];  
            data.chls=[data.chls, chls]; 
            data.isline=[data.isline, isline]; 
            data.n=data.n+1;
            set(handles.axes1, 'UserData', data);  
            
            % Refresh line and point numbers
            RefreshPlot_Callback(hObject, callbackdata, handles);            
            
            % Empty the list of deleted curves
            Reset(hObject, callbackdata, handles, argin);    
        case 'LinePT'
            hedit2=argin{2}; 
            hedit1=argin{4};     
            hedit3=argin{5};     
            
            % Get the point and the length of the line
            pnt1=str2num(get(hedit1, 'String'));
            pnt1=GetCoords(handles, pnt1);
            Len=str2num(get(hedit3, 'String'));
            
            % Get the normal line
            data=get(handles.axes1,'UserData'); 
            lnum=str2num(get(hedit2, 'String')); 
            t=data.LNumber(lnum);
            crv=data.crvs(t);
            
            % Get the normal at the point
            pnt1=[pnt1, 0]';
            dr=neartangent(crv, pnt1);
            dt=[-dr(2); dr(1); dr(3)];
            pnt2=pnt1+Len*dt;
            crvs=nrlline(pnt1', pnt2');
            chls = nrlcrvplot(crvs, 50);   

            % Set the axes limite
            XLim=get(handles.axes1,'XLim');
            YLim=get(handles.axes1,'YLim');
            set(handles.editXLim,'String', num2str(XLim));
            set(handles.editYLim,'String', num2str(YLim));

            % Select curves for further operation
            set(chls, 'ButtonDownFcn', {@Selected, handles, {'Selected'}});

            % Judge whether the entity is a line or a point
            n=length(chls);
            isline=ones(1,n);

            % Set Tag to line handls
            m=length(data.crvs);
            for i=1:n
                set(chls(i), 'Tag', num2str(m+i));
            end

            % Save data to 'UserData' of handles.axes1
            data.crvs=[data.crvs, crvs];  
            data.chls=[data.chls, chls]; 
            data.isline=[data.isline, isline]; 
            data.n=data.n+1;
            set(handles.axes1, 'UserData', data);  
            
            % Refresh line and point numbers
            RefreshPlot_Callback(hObject, callbackdata, handles);            
            
            % Empty the list of deleted curves
            Reset(hObject, callbackdata, handles, argin); 
        case 'RePoints'
            hedit1=argin{2};
            
            % Get the choice of removing points
            dpi=str2num(get(hedit1, 'String'));
            
            % Get the points
            data=get(handles.axes1,'UserData'); 
            crvs=data.crvs;
            n=numel(crvs);
            pnts=[]; t=1; T=[];
            for i=1:n
                if strcmp(crvs(i).form, 'POINTS')
                    pnts=[pnts, crvs(i).coefs'];
                    T(t)=i;
                    t=t+1;
                end
            end
            
            % Find and remove dumplicated points 
            [~, uu]=RemDuplicate(pnts');
            if dpi
                T=T(~uu);
            end
            if ~isempty(T)
                delete(data.chls(T)); 
                data.dcrvs=[data.dcrvs, data.crvs(T)]; 
                data.disline=[data.disline, data.isline(T)]; 
                data.crvs(T)=[]; 
                data.chls(T)=[]; 
                data.isline(T)=[]; 
                set(handles.axes1, 'UserData', data); 
            end
            
            % Empty the list of deleted curves
            Reset(hObject, callbackdata, handles, argin);
            
            % Refresh line and point numbers
            RefreshPlot_Callback(hObject, callbackdata, handles);
            
            % Assign function handle to point numbers and curve handles
            data=get(handles.axes1,'UserData'); 
            set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, argin});
            set(handles.axes1,'UserData', data);       
    end
    
    
% Remove dumplicated points
function ReDpPoints(hObject, callbackdata, handles, argin)
% Get the points
data=get(handles.axes1,'UserData'); 
crvs=data.crvs;
n=numel(crvs);
pnts=[]; t=1; T=[];
for i=1:n
    if strcmp(crvs(i).form, 'POINTS')
        pnts=[pnts, crvs(i).coefs'];
        T(t)=i;
        t=t+1;
    end
end

% Find and remove dumplicated points 
[~, uu]=RemDuplicate(pnts');
T=T(~uu);
if ~isempty(T)
    delete(data.chls(T)); 
    data.dcrvs=[data.dcrvs, data.crvs(T)]; 
    data.disline=[data.disline, data.isline(T)]; 
    data.crvs(T)=[]; 
    data.chls(T)=[]; 
    data.isline(T)=[]; 
    set(handles.axes1, 'UserData', data); 
end

% Empty the list of deleted curves
Reset(hObject, callbackdata, handles, argin);

% Refresh line and point numbers
RefreshPlot_Callback(hObject, callbackdata, handles);

% Assign function handle to point numbers and curve handles
data=get(handles.axes1,'UserData'); 
set(data.hPNum, 'ButtonDownFcn', {@Selected, handles, argin});
set(handles.axes1,'UserData', data);      
    
    
% Delete function for uihandles
function DeleteFcn(hObject, callbackdata, handles, argin) 
    switch argin{1}
        case {'Selected', 'Mouse'}
            data=get(handles.axes1,'UserData'); 
            set(data.chls, 'ButtonDownFcn', {@Selected, handles, {'Selected'}});     
            set(gcf, 'WindowButtonDownFcn', {@SetCurrentCoords, handles, {'none'}});
    end

% --- Executes on button press in RePointsButton.
function RePointsButton_Callback(hObject, eventdata, handles)
% hObject    handle to RePointsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
RemovePoints_Callback(hObject, eventdata, handles);


% --- Executes on button press in MeasureButton.
function MeasureButton_Callback(hObject, eventdata, handles)
% hObject    handle to MeasureButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Measure_Callback(hObject, eventdata, handles);


% --- Executes on button press in Recover.
function Recover_Callback(hObject, eventdata, handles)
% hObject    handle to Recover (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.axes1,'UserData'); 
nd=numel(data.dcrvs);
if nd>0
    crvs=data.dcrvs(nd);
    if data.disline(nd)
        chls= nrbcrvplot(crvs,50); 
    else
        chls=pntplot(pnt, 0.5, '.');
    end
    
    % Set the axes limite
    XLim=get(handles.axes1,'XLim');
    YLim=get(handles.axes1,'YLim');
    set(handles.editXLim,'String', num2str(XLim));
    set(handles.editYLim,'String', num2str(YLim));
    
    % Select curves for further operation
    set(chls, 'ButtonDownFcn', {@Selected, handles, {'Selected'}});
    
    % Judge whether the entity is a line or a point
    data=get(handles.axes1, 'UserData');  
    n=length(chls);
    if ~data.disline(nd)
        isline=0;
    else
        isline=ones(1,n);
    end
    
    % Set Tag to line handls
    m=length(data.crvs);
    for i=1:n
        set(chls(i), 'Tag', num2str(m+i));
    end
        
    % Delete the curve from 'UserData'
    data.dcrvs(nd)=[];
    data.disline(nd)=[];    
    
    % Save data to 'UserData' of handles.axes1
    data.crvs=[data.crvs, crvs];  
    data.chls=[data.chls, chls]; 
    data.isline=[data.isline, isline]; 
    data.n=data.n+1;    
    set(handles.axes1, 'UserData', data);  
end
