function varargout = ms_calc(varargin)
% MS_CALC MATLAB code for ms_calc.fig
%  MS_CALC, by itself, creates a new MS_CALC or raises the existing
%  singleton*.
%
%  H = MS_CALC returns the handle to a new MS_CALC or the handle to
%  the existing singleton*.
%
%  MS_CALC('CALLBACK',hObject,eventData,handles,...) calls the local
%  function named CALLBACK in MS_CALC.M with the given input arguments.
%
%  MS_CALC('Property','Value',...) creates a new MS_CALC or raises the
%  existing singleton*.  Starting from the left, property value pairs are
%  applied to the GUI before ms_calc_OpeningFcn gets called.  An
%  unrecognized property name or invalid value makes property application
%  stop.  All inputs are passed to ms_calc_OpeningFcn via varargin.
%
%  *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%  instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ms_calc

% Last Modified by GUIDE v2.5 27-Mar-2017 20:47:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',   mfilename, ...
'gui_Singleton',  gui_Singleton, ...
'gui_OpeningFcn', @ms_calc_OpeningFcn, ...
'gui_OutputFcn',  @ms_calc_OutputFcn, ...
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


% --- Executes just before ms_calc is made visible.
function ms_calc_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObjecthandle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesstructure with handles and user data (see GUIDATA)
% varargin   command line arguments to ms_calc (see VARARGIN)

% Choose default command line output for ms_calc
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ms_calc wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ms_calc_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObjecthandle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesstructure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObjecthandle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesstructure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObjecthandle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesempty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%   See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function edit_episilon_Callback(hObject, eventdata, handles)
% hObjecthandle to edit_episilon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesstructure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_episilon as text
%str2double(get(hObject,'String')) returns contents of edit_episilon as a double


% --- Executes during object creation, after setting all properties.
function edit_episilon_CreateFcn(hObject, eventdata, handles)
% hObjecthandle to edit_episilon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesempty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%   See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function edit_h_Callback(hObject, eventdata, handles)
% hObjecthandle to edit_h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesstructure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_h as text
%str2double(get(hObject,'String')) returns contents of edit_h as a double


% --- Executes during object creation, after setting all properties.
function edit_h_CreateFcn(hObject, eventdata, handles)
% hObjecthandle to edit_h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesempty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%   See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function edit_cond_Callback(hObject, eventdata, handles)
% hObjecthandle to edit_cond (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesstructure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cond as text
%str2double(get(hObject,'String')) returns contents of edit_cond as a double


% --- Executes during object creation, after setting all properties.
function edit_cond_CreateFcn(hObject, eventdata, handles)
% hObjecthandle to edit_cond (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesempty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%   See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in pb_synthesis.
function pb_synthesis_Callback(hObject, eventdata, handles)
% hObjecthandle to pb_synthesis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesstructure with handles and user data (see GUIDATA)
h = str2num(get(handles.edit_h,'string')); %substrate thickness
%w = float(lineEdit_ms_W.text())/1.0e3 % conductor width
%l = float(lineEdit_ms_L.text())/1.0e3 % conductor length
c = 2.99792458e8;
er =  str2num(get(handles.edit_episilon,'string')) ;% relative dielectric constant
mur =  str2num(get(handles.edit_mu,'string')); % relative permeability
fc =  str2num(get(handles.edit_fc,'string'))*1e6; % frequency in Hz
rough = str2num(get(handles.edit_rought,'string'))/1.0e3;
tand = str2num(get(handles.edit_tand,'string')); % loss tangent of the dielectric
cond = str2num(get(handles.edit_cond,'string'));
mu = mur*4*pi*1.0e-7;
Zc = str2num(get(handles.edit_Zc,'string'));
t = str2num(get(handles.edit_T,'string'))/1.0e3;
el = str2num(get(handles.edit_eL,'string'))/1.0e3;
lx = 1000*25.4e-6;
wmin = 0.01*25.4e-6;wmax = 499.0*25.4e-6;
%impedance convergence tolerance (ohms)
abstol = 1e-6;
reltol = 0.1e-6;
maxiters = 50;
A = ((er - 1)/(er + 1)) * (0.226 + 0.121/er) + (pi/377.0)*sqrt(2*(er+1))*Zc;
w_h = 4/(0.5*exp(A) - exp(-A));
if w_h > 2.0 
B = pi*377.0/(2*Zc*sqrt(er));
w_h = (2/pi)*(B - 1 - log(2*B - 1) + ((er-1)/(2*er))*(log(B-1) + 0.293 - 0.517/er));
end
wx = h * w_h;
if wx >= wmax
wx = 0.95*wmax;
end
if wx <= wmin
wx = wmin;
end
wold = 1.01*wx;
Zold = z_calc(wold,h,lx,t,fc,er);
if Zold < Zc
wmax = wold;
else
wmin = wold;
end
iters = 0;
done = 0;
while done == 0
iters = iters + 1;
Z0 = z_calc(wx,h,lx,t,fc,er);
if Z0 < Zc
wmax = wx;
else
wmin = wx;
end
if abs(Z0-Zc) < abstol
done = 1;
elseif abs(wx-wold) < reltol
done = 1;
elseif iters >=  maxiters 
done = 1;
else
dzdw = (Z0 -Zold)/(wx-wold);
wold = wx;
Zold = Z0;
wx = wx -(Z0-Zc)/dzdw;
if (wx > wmax) |(wx < wmin)
wx = (wmin+ wmax)/2.0;
end
end
end
er_eff = er_eff_calc(wx,h,t,fc,er);
v = c/sqrt(er_eff);
l = el/360.0*v/fc;
%lambda0= c/fc; %m
msa_w = c*1000.0/2.0/fc*sqrt(2.0/(er+1.0)); % patch width  m
set(handles.edit_W,'string',num2str(msa_w));
end
 


% --- Executes on button press in pb_analysis.
function pb_analysis_Callback(hObject, eventdata, handles)
% hObjecthandle to pb_analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesstructure with handles and user data (see GUIDATA)



function edit_tand_Callback(hObject, eventdata, handles)
% hObjecthandle to edit_tand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesstructure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tand as text
%str2double(get(hObject,'String')) returns contents of edit_tand as a double


% --- Executes during object creation, after setting all properties.
function edit_tand_CreateFcn(hObject, eventdata, handles)
% hObjecthandle to edit_tand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesempty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%   See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function edit_T_Callback(hObject, eventdata, handles)
% hObjecthandle to edit_T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesstructure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_T as text
%str2double(get(hObject,'String')) returns contents of edit_T as a double


% --- Executes during object creation, after setting all properties.
function edit_T_CreateFcn(hObject, eventdata, handles)
% hObjecthandle to edit_T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesempty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%   See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function edit_rought_Callback(hObject, eventdata, handles)
% hObjecthandle to edit_rought (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesstructure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rought as text
%str2double(get(hObject,'String')) returns contents of edit_rought as a double

% --- Executes during object creation, after setting all properties.
function  edit_rought_CreateFcn(hObject, eventdata, handles)
% hObjecthandle to edit_rought (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesempty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%   See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end

function  edit_Zc_Callback(hObject, eventdata, handles)

% hObjecthandle to edit_Zc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesstructure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Zc as text
%str2double(get(hObject,'String')) returns contents of edit_Zc as a double

% --- Executes during object creation, after setting all properties.
function edit_Zc_CreateFcn(hObject, eventdata, handles)
% hObjecthandle to edit_Zc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesempty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%   See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end


function edit_mu_Callback(hObject, eventdata, handles)
% hObjecthandle to edit_mu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesstructure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mu as text
%str2double(get(hObject,'String')) returns contents of edit_mu as a double


% --- Executes during object creation, after setting all properties.
function edit_mu_CreateFcn(hObject, eventdata, handles)
% hObjecthandle to edit_mu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesempty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%   See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function edit_L_Callback(hObject, eventdata, handles)
% hObjecthandle to edit_L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesstructure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_L as text
%str2double(get(hObject,'String')) returns contents of edit_L as a double


% --- Executes during object creation, after setting all properties.
function edit_L_CreateFcn(hObject, eventdata, handles)
% hObjecthandle to edit_L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesempty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%   See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function edit_W_Callback(hObject, eventdata, handles)
% hObjecthandle to edit_W (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesstructure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_W as text
%str2double(get(hObject,'String')) returns contents of edit_W as a double


% --- Executes during object creation, after setting all properties.
function edit_W_CreateFcn(hObject, eventdata, handles)
% hObjecthandle to edit_W (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesempty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%   See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function edit_fc_Callback(hObject, eventdata, handles)
% hObjecthandle to edit_fc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesstructure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fc as text
%str2double(get(hObject,'String')) returns contents of edit_fc as a double


% --- Executes during object creation, after setting all properties.
function edit_fc_CreateFcn(hObject, eventdata, handles)
% hObjecthandle to edit_fc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesempty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%   See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function edit_eL_Callback(hObject, eventdata, handles)
% hObjecthandle to edit_eL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesstructure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_eL as text
%str2double(get(hObject,'String')) returns contents of edit_eL as a double


% --- Executes during object creation, after setting all properties.
function edit_eL_CreateFcn(hObject, eventdata, handles)
% hObjecthandle to edit_eL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesempty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%   See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function edit_beta_Callback(hObject, eventdata, handles)
% hObjecthandle to edit_beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesstructure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_beta as text
%str2double(get(hObject,'String')) returns contents of edit_beta as a double


% --- Executes during object creation, after setting all properties.
function edit_beta_CreateFcn(hObject, eventdata, handles)
% hObjecthandle to edit_beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesempty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%   See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function edit_episilon_e_Callback(hObject, eventdata, handles)
% hObjecthandle to edit_episilon_e (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesstructure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_episilon_e as text
%str2double(get(hObject,'String')) returns contents of edit_episilon_e as a double


% --- Executes during object creation, after setting all properties.
function edit_episilon_e_CreateFcn(hObject, eventdata, handles)
% hObjecthandle to edit_episilon_e (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesempty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%   See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end



function edit_loss_Callback(hObject, eventdata, handles)
% hObjecthandle to edit_loss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesstructure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_loss as text
%str2double(get(hObject,'String')) returns contents of edit_loss as a double


% --- Executes during object creation, after setting all properties.
function edit_loss_CreateFcn(hObject, eventdata, handles)
% hObjecthandle to edit_loss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesempty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%   See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObjecthandle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlesstructure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1

function Zc = z_calc(w,h,l,t,fc,er)
%fn = fc*h/1e6
U = w/h ;% ratio of trace width to substrate thickness 
if t > 0
T = t/h; %ratio of conductor thickness to substrate thickness
%(T/PI)*log(1.0 + 4.0*exp(1.0)/(T*power(coth(sqrt(6.517*u)),2.0)))
U1 = U +(T*log(1.0+4.0*e/T/power(1.0/tanh(sqrt(6.517*U)),2.0)))/pi; % from Hammerstad and Jensen
%   0.5*(1.0 + 1.0/cosh(sqrt(er-1.0)))*deltau1
Ur = U +(U1-U)*(1.0+1.0/(cosh(sqrt(er-1))))/2.0; % from Hammerstad and Jensen
else
U1 = U ; 
Ur = U;
end
Y = ee_HandJ(Ur,er);
Z0 =z0_HandJ(Ur)/sqrt(Y);
%ereff0 = Y*power(Z01_U1/Z01_Ur,2)
ereff0 = Y*power(z0_HandJ(U1)/z0_HandJ(Ur),2.0);
fn = fc*h/1e7;
P1 = 0.27488 + (0.6315 + (0.525 / (power((1 + 0.157*fn),20))) )*U - 0.065683*exp(-8.7513*U);
P2 = 0.33622*(1 - exp(-0.03442*er));
P3 = 0.0363*exp(-4.6*U)*(1 - exp(-power((fn / 3.87),4.97)));
P4 = 1 + 2.751*( 1 -  exp(-power((er/15.916),8)));
P = P1*P2*power(((0.1844 + P3*P4)*10*fn),1.5763);
ereff = (er*P+ereff0)/(1+P); % equavlent ralative dielectric constant
%ms_Zc = Z0*power(ereff0/ms_ereff,0.5)*(ms_ereff-1)/(ereff0-1)
fn = fc*h/1e6;
R1 = 0.03891*(power(er,1.4));
R2 = 0.267*(power(U,7.0));
R3 = 4.766*exp(-3.228*(power(U,0.641)));
R4 = 0.016 + power((0.0514*er),4.524);
R5 = power((fn/28.843),12.0);
R6 = 22.20*(power(U,1.92));
R7 = 1.206 - 0.3144*exp(-R1)*(1 - exp(-R2));
R8 = 1.0 + 1.275*(1.0 -  exp(-0.004625*R3*power(er,1.674)*power(fn/18.365,2.745)));
R9 = (5.086*R4*R5/(0.3838 + 0.386*R4))*(exp(-R6)/(1 + 1.2992*R5));
R9 = R9 * (power((er-1),6))/(1 + 10*power((er-1),6));
R10 = 0.00044*(power(er,2.136)) + 0.0184;
R11 = (power((fn/19.47),6))/(1 + 0.0962*(power((fn/19.47),6)));
R12 = 1 / (1 + 0.00245*U*U);
R13 = 0.9408*(power( ereff,R8)) - 0.9603;
R14 = (0.9408 - R9)*(power(ereff0,R8))-0.9603;
R15 = 0.707*R10*(power((fn/12.3),1.097));
R16 = 1 + 0.0503*er*er*R11*(1 - exp(-(power((U/15),6))));
R17 = R7*(1 - 1.1241*(R12/R16)*exp(-0.026*(power(fn,1.15656))-R15));
Zc = Z0*(power((R13/R14),R17));% characteristic impedance
end

function er_eff =er_eff_calc(w,h,t,fc,er)
%fn = fc*h/1e6
U = w/h; % ratio of trace width to substrate thickness 
if t > 0
T = t/h; %ratio of conductor thickness to substrate thickness
%(T/PI)*log(1.0 + 4.0*exp(1.0)/(T*power(coth(sqrt(6.517*u)),2.0)))
U1 = U +(T*log(1.0+4.0*e/T/power(1.0/tanh(sqrt(6.517*U)),2.0)))/pi; % from Hammerstad and Jensen
%   0.5*(1.0 + 1.0/cosh(sqrt(er-1.0)))*deltau1
Ur = U +(U1-U)*(1.0+1.0/(cosh(sqrt(er-1))))/2.0; % from Hammerstad and Jensen
else
U1 = U;  
Ur = U;
end
Y = ee_HandJ(Ur,er);
Z0 =z0_HandJ(Ur)/sqrt(Y);
%ereff0 = Y*power(Z01_U1/Z01_Ur,2)
ereff0 = Y*power(z0_HandJ(U1)/z0_HandJ(Ur),2.0);
fn = fc*h/1e7;
P1 = 0.27488 + (0.6315 + (0.525 / (power((1 + 0.157*fn),20))) )*U - 0.065683*exp(-8.7513*U);
P2 = 0.33622*(1 - exp(-0.03442*er));
P3 = 0.0363*exp(-4.6*U)*(1 - exp(-power((fn / 3.87),4.97)));
P4 = 1 + 2.751*( 1 -  exp(-power((er/15.916),8)));
P = P1*P2*power(((0.1844 + P3*P4)*10*fn),1.5763);
er_eff = (er*P+ereff0)/(1+P);% equavlent ralative dielectric constant
%return er_eff
end

function z01 = z0_HandJ(u)
LIGHTSPEED = 2.99792458e8;
FREESPACEZ0 = 4.0*pi*1.0e-7*LIGHTSPEED;
%from Hammerstad and Jensen.  'u' is the normalized width
F = 6.0 + (2.0*pi - 6.0)*exp(-1*power((30.666/u),0.7528));
%from Hammerstad and Jensen
z01 = (FREESPACEZ0/(2*pi))*log(F/u + sqrt(1.0 + power((2/u),2.0)));
%return z01
end

function Y = ee_HandJ(u,er)
%Au = 1.0 + (1.0/49.0)*log(power(Ur,4.0) + power(Ur/52.0,2.0)/(power(Ur,4.0)+0.432))+ (1.0/18.7)*log(1.0+ power(Ur/18.1,3.0))
A = 1.0 + (1.0/49.0)*log((power(u,4.0) + power((u/52.0),2.0))/(power(u,4.0) + 0.432))...
+ (1.0/18.7)*log(1.0 + power((u/18.1),3.0));
%Ber = 0.564*power((er-0.9)/(er+3),0.053)
B = 0.564*power(((er-0.9)/(er+3.0)),0.053);
%Y = (er+1.0)/2.0+(er-1.0)/2.0*power(1+10.0/Ur,-(A*B))
Y= (er+1.0)/2.0 + ((er-1.0)/2.0)*power((1.0 + 10.0/u),(-A*B));
end


