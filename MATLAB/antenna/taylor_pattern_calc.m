function varargout = taylor_pattern_calc(varargin)
% TAYLOR_PATTERN_CALC MATLAB code for taylor_pattern_calc.fig
%      TAYLOR_PATTERN_CALC, by itself, creates a new TAYLOR_PATTERN_CALC or raises the existing
%      singleton*.
%
%      H = TAYLOR_PATTERN_CALC returns the handle to a new TAYLOR_PATTERN_CALC or the handle to
%      the existing singleton*.
%
%      TAYLOR_PATTERN_CALC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TAYLOR_PATTERN_CALC.M with the given input arguments.
%
%      TAYLOR_PATTERN_CALC('Property','Value',...) creates a new TAYLOR_PATTERN_CALC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before taylor_pattern_calc_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to taylor_pattern_calc_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help taylor_pattern_calc

% Last Modified by GUIDE v2.5 25-Apr-2017 15:43:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @taylor_pattern_calc_OpeningFcn, ...
    'gui_OutputFcn',  @taylor_pattern_calc_OutputFcn, ...
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

% --- Executes just before taylor_pattern_calc is made visible.
function taylor_pattern_calc_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to taylor_pattern_calc (see VARARGIN)

% Choose default command line output for taylor_pattern_calc
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using taylor_pattern_calc.
if strcmp(get(hObject,'Visible'),'off')
    plot(rand(5));
end

% UIWAIT makes taylor_pattern_calc wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = taylor_pattern_calc_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
cla;
global amp;
global N;
global SLL0;
popup_sel_index = get(handles.popupmenu1, 'Value');
N = str2num(get(handles.edit_elem,'String'));
SLL0 = str2num(get(handles.edit_SLL,'String'));
switch popup_sel_index
    case 1
        %         plot(rand(5));
        
        B1= 0.1;
        B2 = 2;
        SLL1= SLL_CALC(B1);
        SLL2= SLL_CALC(B2);
        while(abs(SLL1+SLL0) > 0.001 || abs(SLL2+SLL0)> 0.001)
            B =(B1+B2)/2;
            SLL = SLL_CALC(B);
            if SLL + SLL0 > 0
                B2 =B;
            else
                B1=B;
            end
            SLL1= SLL_CALC(B1);
            SLL2= SLL_CALC(B2);
        end
        B =(B1+B2)/2;
        SLL= SLL_CALC(B);
        
        lambda = 1;
        d = lambda/2;
        L = d*(N-1);
        i = zeros(1,N);
        
        for m = 1:1:N
            i(m) = besselj(0,1i*pi*B*sqrt(1-(2*(m-(N+1)/2)*d/L)^2));
        end
        amp = i/max(i)
        plot(amp);
         grid on;
        xlim([1 N]);
        ylim([min(amp) 1]);
        str = strcat('N=', num2str(N),', SLL = ',num2str(SLL0) ,'dB   Taylor 单变量综合');
        title(str);
        
    case 2
        R0 = 10^(-SLL0/20);
        x = zeros(1,N-1);
        u = zeros(1,N-1);
        s = zeros(1,N);
        I = zeros(1,N);
        x0 = 1/2*(power(R0+sqrt(R0^2-1),1/(N-1))+power(R0-sqrt(R0^2-1),1/(N-1)));
        if mod(N,2) == 1
            K = round((N-1)/2);
            
            for n =2:1:K + 1
                for p = n:1:K+1
                    I(n) = I(n)+(-1)^(K-p+1)*K*factorial(p+K-2)/factorial(p-n)/factorial(p+n-2)/factorial(K-p+1)*power(x0,2*(p-1))   ;
                end
                
            end
            for p = 1:1:K+1
                I(1) = I(1)+(-1)^(K-p+1)*K*factorial(p+K-2)/factorial(p-1)/factorial(p-1)/factorial(K-p+1)*power(x0,2*(p-1))   ;
            end
            
            
            %x(1) = 1;
            %x(N) = 1;
            for i = 1:1:K+1
                x(i) = I(K+2-i)/I(K+1);
                x(N+1-i) = x(i);
            end
            
        else
            K = round(N/2);
            
            for n =1:1:K
                for p = n:1:K
                    I(n) = I(n)+(-1)^(K-p)*(2*K-1)*factorial(p+K-2)/2/factorial(p-n)/factorial(p+n-1)/factorial(K-p)*power(x0,2*p-1);
                end
                I(N+1-n) = I(n);
            end
            
            for i = 1:1:K
                x(i) = I(K+1-i)/I(K);
                x(N+1-i) = x(i);
            end
        end
        amp = x/max(x);
        plot(amp);
         grid on;
        xlim([1 N]);
        ylim([min(amp) 1]);
        str = strcat('N=', num2str(N),', SLL = ',num2str(SLL0) ,'dB by Dolph-Chebyshev Synthesis');
        title(str);
    case 3
        amp = ones(1,N);
         plot(amp);
          grid on;
        xlim([1 N]);
        ylim([0.5 1.5]); 
        str = strcat('N=', num2str(N),' 均匀分布');
        title(str);
%     case 4
%         plot(membrane);
%     case 5
%         surf(peaks);
end
set(handles.togglebutton,'String','Pattern');


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
    ['Close ' get(handles.figure1,'Name') '...'],...
    'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'Taylor 单变量', 'Dolph-Chebyshev)', '均匀分布'});


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_elem_Callback(hObject, eventdata, handles)
% hObject    handle to edit_elem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_elem as text
%        str2double(get(hObject,'String')) returns contents of edit_elem as a double


% --- Executes during object creation, after setting all properties.
function edit_elem_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_elem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_SLL_Callback(hObject, eventdata, handles)
% hObject    handle to edit_SLL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_SLL as text
%        str2double(get(hObject,'String')) returns contents of edit_SLL as a double


% --- Executes during object creation, after setting all properties.
function edit_SLL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SLL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in togglebutton.
function togglebutton_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton
global amp;
global N;
global SLL0;
amp
axes(handles.axes1);
cla;
popup_sel_index = get(handles.popupmenu1, 'Value');
if strcmp(get(handles.togglebutton,'String'),'Pattern')
    set(handles.togglebutton,'String','Amp');
    
    x = amp;
    M = 2000; % theta sample points
    E0 = zeros(1,M);%
    H0 = zeros(1,M);%
    U0= zeros(1,M);%
    f0 = 1e9;
    lambda = 3e8/f0;
    d = lambda/2;% elements spacing
    k =2*pi/lambda;% wave constant
    imp = 120*pi;%wave impdance
    I0= 1;
    r = 100;
    l = 0.01;
    
    af = zeros(1,M);% array factor buffer
    
    theta = linspace(0,pi,M)+pi/10e10;
    phi = linspace(0,2*pi,M);
    af = zeros(1,M);
    a = 0;
    b = 0;
    for j = 1:1:N
        a= a+x(j);
        b = b + power(x(j),2);
        af = af+ exp(1i*(j-1)*k*d*cos(theta))*x(j);
    end
    af = abs(af/N );
    E0 = af *imp;
    Wav = real(E0 .* E0 )/2/imp;
    Prad = imp*pi*3*(I0*l/lambda)^2;
    theta = theta/pi*180;
    af= 10*log10(af/max(af));
    U= real(E0 .*E0)*r^2/2/imp;
    U= 10*log10(U/max(U));
    
    error = 0.01;
    af=U;
    pos1_3dB = [];
    pos_max = find(max(af)==af);
    while(isempty(pos1_3dB))
        pos1_3dB = find(abs(((af(1:pos_max)-af(pos_max)))+3) < error);
        error = error + 0.005;
    end
    error = 0.01;
    pos2_3dB = [];
    while(isempty(pos2_3dB))
        pos2_3dB = find(abs(((af(pos_max:end)-af(pos_max)))+3) < error);
        error = error + 0.005;
    end
    BeamWidth= (theta(pos2_3dB(1)+pos_max)-theta(pos1_3dB(end)));
    D = 10*log10(a^2/b);
    tau = a^2/b/N;
    
    plot(theta,U);
    grid on;
    hold on;
    axis([0 180 -120 0]);
    %             str = strcat('N=', num2str(N),', SLL = ',num2str(SLL0) ,'dB by Taylor Synthesis');
    bw = strcat('main lobe beamwidth=',num2str(BeamWidth),'degree, gain =  ',num2str(D),'dB')
    %             title(str)
    text(60,-10,bw,'fontsize',10);
    ylim([-60 0])
    
else
    set(handles.togglebutton,'String','Pattern');
    plot(amp);
     grid on;
    xlim([1 N]);
    if popup_sel_index ~= 3
    ylim([min(amp) 1]);
    else 
        ylim([0.5 1.5]);
    end
end

switch popup_sel_index
    case 1
        str = strcat('N=', num2str(N),', SLL = ',num2str(SLL0) ,'dB by Taylor Synthesis');
        title(str);
         
    case 2
        str = strcat('N=', num2str(N),', SLL = ',num2str(SLL0) ,'dB by Dolph-Chebyshev  Synthesis');
        title(str);
        
    case 3
        str = strcat('N=', num2str(N),' 均匀分布');
        title(str);
        
end