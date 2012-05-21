function varargout = gui_inp(varargin)
% GUI_INP M-file for gui_inp.fig
%      GUI_INP, by itself, creates a new GUI_INP or raises the existing
%      singleton*.
%
%      H = GUI_INP returns the handle to a new GUI_INP or the handle to
%      the existing singleton*.
%
%      GUI_INP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_INP.M with the given input arguments.
%
%      GUI_INP('Property','Value',...) creates a new GUI_INP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_inp_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_inp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_inp

% Last Modified by GUIDE v2.5 12-Aug-2010 15:34:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @gui_inp_OpeningFcn, ...
    'gui_OutputFcn',  @gui_inp_OutputFcn, ...
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

function gui_inp_OpeningFcn(hObject, eventdata, handles, varargin)
global file_dir file_name inp_co_def plan_lmts txt looplimits

file_dir = pwd;                 %current working dir
b='\config_file.txt';           %default file name
file_name=strcat(file_dir,b);   %overall path name
set(handles.saveloc_edit,'String',file_name)    %display in the path text box GUI

s=xlsread('defaults','INP_distrb')';    %INP distributions from defaults excel sheet
inp_co_def(1:16)=s(2,:);                %second column of the INP_distrb worksheet
T_sum=sum(inp_co_def);                  %Total sum
inp_co_def(:)=inp_co_def(:)./T_sum;     %normalizing values

[plan_lmts,txt]=xlsread('defaults','plan');         %read all possible plans combinations via defaults.xls
set(handles.plan_set_menu,'String',txt(2:end,1));   %setting the 'freq plan' strings in GUI
set(handles.remarks,'String',txt(2,2));             %setting 'remark' strings of the freq plan in GUI

plan_chosen=1;                              %when initialized, first plan is chosen as default
looplimits=zeros(1,23);                        %initializing

looplimits(20)=plan_lmts(plan_chosen,1);    %plan chosen adsl/vdsl -opening defaults
looplimits(22)=plan_lmts(plan_chosen,2);    %downstream/upstream

handles.output = hObject;
guidata(hObject, handles);

function plan_set_menu_Callback(hObject, eventdata, handles)

global looplimits plan_lmts txt

plan_chosen=get(handles.plan_set_menu,'Value');     %plan chosen from the pop menu
set(handles.remarks,'String',txt(plan_chosen+1,2)); %remarks display -GUI

if plan_lmts(plan_chosen,1)==1                      %vdsl case
    set(handles.Lpmin_edit,'String','32')           %corresponding limits
    set(handles.Lpmax_edit,'String','30000')
    set(handles.Qmin_edit,'Enable','On')
    set(handles.Qmax_edit,'Enable','On')
    set(handles.Qmax_edit,'String','8')
else set(handles.Lpmin_edit,'String','16')          %adsl case
    set(handles.Lpmax_edit,'String','7000')         %corresponding limits
    set(handles.Qmax_edit,'String','1')
    set(handles.Qmin_edit,'Enable','Off')
    set(handles.Qmax_edit,'Enable','Off')
end

looplimits(14)=plan_lmts(plan_chosen,7);    %Dmax
looplimits(15)=plan_lmts(plan_chosen,5);    %inv smin
looplimits(16)=plan_lmts(plan_chosen,6);    %smax
looplimits(17)=plan_lmts(plan_chosen,3);    %max_intv
looplimits(18)=plan_lmts(plan_chosen,4);    %Fs

looplimits(19)=str2num(get(handles.msg_edit,'String')); %msg_min

looplimits(20)=plan_lmts(plan_chosen,1);                %adsl/vdsl?

checkboxStatus = get(handles.force_check,'Value');      %fast path case
if (checkboxStatus) 
    looplimits(21)=1;           %fast path
else looplimits(21)=2;          %full path
end

looplimits(22)=plan_lmts(plan_chosen,2);                %downstream/upstream

checkboxStatus = get(handles.framing_check,'Value');    %framing
if (checkboxStatus) 
    looplimits(23)=1;           %framing
else looplimits(23)=2;          %no framing
end

function start_button_Callback(hObject, eventdata, handles)

global looplimits file_name plan_lmts   %global variables
global D_uv Gp_uv Tp_uv N_uv Mp_uv

looplimits(1:23)=0;                     %Initializing to zero.

set(handles.dummy,'String','Reading');  %Status bar text
    
plan_chosen=get(handles.plan_set_menu,'Value'); %the plan chosen from the drop box

looplimits(1)=str2num(get(handles.INPmin_edit,'String'));   %limits of INP
looplimits(2)=str2num(get(handles.INPmax_edit,'String'));
looplimits(3)=str2num(get(handles.INPstep_edit,'String'));

looplimits(4)=str2num(get(handles.delaymin_edit,'String')); %limits of delay
looplimits(5)=str2num(get(handles.delaymax_edit,'String'));
looplimits(6)=str2num(get(handles.delaystep_edit,'String'));

looplimits(7)=str2num(get(handles.Lpmin_edit,'String'));    %limits of Lp
looplimits(8)=str2num(get(handles.Lpmax_edit,'String'));
looplimits(9)=str2num(get(handles.Lpstep_edit,'String'));

looplimits(10)=str2num(get(handles.Rmin_edit,'String'));    %limits of R
looplimits(11)=str2num(get(handles.Rmax_edit,'String'));

looplimits(12)=str2num(get(handles.Qmin_edit,'String'));    %limits of Q
looplimits(13)=str2num(get(handles.Qmax_edit,'String'));

looplimits(14)=plan_lmts(plan_chosen,7);    %Dmax
looplimits(15)=plan_lmts(plan_chosen,5);    %Inv smin
looplimits(16)=plan_lmts(plan_chosen,6);    %Smax
looplimits(17)=plan_lmts(plan_chosen,3);    %Max_intv
looplimits(18)=plan_lmts(plan_chosen,4);    %Fs

looplimits(19)=str2num(get(handles.msg_edit,'String'));     %msg_min

looplimits(20)=plan_lmts(plan_chosen,1);                    %plan chosen adsl/vdsl

checkboxStatus = get(handles.force_check,'Value');          %fast path case
if (checkboxStatus) 
    looplimits(21)=1;                                       %fast path case D=1 case
else looplimits(21)=2;                                      %full path 
end

looplimits(22)=plan_lmts(plan_chosen,2);                    %downstream/upstream

checkboxStatus = get(handles.framing_check,'Value');        %framing check
if (checkboxStatus) 
    looplimits(23)=1;       %framing
else looplimits(23)=2;      %no framing
end

D_adsl=[1,2,4,8,16,32,64,96,128,160,192,224,256,288,320,352,384,416,448,480,511];   %all possible adsl D values
Mp_all=[1,2,4,8,16];                                    %all possible Mp values

file_name=get(handles.saveloc_edit,'String');           %constitutes for possible change in file name

set(handles.dummy, 'String','Busy');                    %status update - GUI
set(handles.start_button,'Enable','Inactive');          %inactivate the start button -GUI

checkboxStatus = get(handles.uniquechkbox,'Value');     %obtain unique values?
if (checkboxStatus)
    if looplimits(20)==1        %vdsl case
        dlim=looplimits(14);    %looplimits(14)=D_max
        disp('Output format(.txt): Inp-L-N-s-D-q-R-msg-Tp-Mp-Gp-perb-per-seq-op');
    elseif looplimits(22)==1
        dlim=21;                %adsl,ds case 21 possible D vals
        disp('Output format(.txt): Inp-L-N-s-D-q-R-msg-Tp-Mp-per-seq-or');
    elseif looplimits(22)==2
        dlim=7;                 %adsl,us case 7 possible D vals
        disp('Output format(.txt): Inp-L-N-s-D-q-R-msg-Tp-Mp-per-seq-or');
    end
    
    %occurance flags
    Gp_uv(1:32)=0;  
    Tp_uv(1:64)=0;
    N_uv(1:255)=0;
    Mp_uv(1:5)=0;
    D_uv(1:dlim)=0;
    
    switch get(handles.unique_popmenu,'Value')  %unique value for which parameter - drop down menu?
        case 1                                          %all unique valus of D
            unique_distb('D');                          %follow INP distrb  
            disp('Done D')              
            if min(D_uv(1:dlim))~=0                     %check if all D met
                disp('Unique D met')
            elseif looplimits(20)==1                    %vdsl 
                D_not_met=find(~D_uv(1:dlim))           %print those not met
            else D_not_met=D_adsl(find(~D_uv(1:dlim)))  %adsl, print those not met.
            end
        case 2                                          %all unique values of N
            unique_distb('N');                          %INP distrb
            disp('Done N')
            if min(N_uv(32:255))~=0                     %check if all N met
                disp('Unique N met')
            else N_not_met = find(~N_uv(32:255))+31     %print those not met
            end
        case 3                                          %all unique Tp values
            unique_distb('Tp');                         %follow INP distrb
            disp('Done Tp')
            if min(Tp_uv(1:64))~=0                      %check if all Tp met
                disp('Unique Tp met')
            else Tp_not_met=find(~Tp_uv(1:64))          %print those Tp not met
            end
        case 4                                          %all unique values of Gp
            unique_distb('Gp');                         %follow distb of INP    
            disp('Done Gp')
            if min(Gp_uv(1:32))~=0                      %check if all Gp met
                disp('Unique Gp met')
            elseif looplimits(20)==1                    %vsdl
                Gp_not_met=find(~Gp_uv(1:32))           %print those not met
            else disp('No Gp for ADSL')                 %adsl, no Gp
            end
        case 5                                          %all unique vals of Mp    
            unique_distb('Mp');                         %follow distrb of INP
            disp('Done Mp')
            if min(Mp_uv(1:5))~=0                       %check if all Mp met
                disp('Unique Mp met')
            else Mp_not_met=Mp_all(find(~Mp_uv(1:5)))   %print those not met
            end
        case 6                                          %Whole unique cases
            unique_distb('D');                          %D
            disp('Done D')
            if min(D_uv(1:dlim))~=0
                disp('Unique D met')
            elseif looplimits(20)==1 %vdsl
                D_not_met=find(~D_uv(1:dlim))
            else D_not_met=D_adsl(find(~D_uv(1:dlim)))
            end
            unique_distb('N');                          %Nfec
            disp('Done N')
            if min(N_uv(32:255))~=0
                disp('Unique N met')
            else N_not_met = find(~N_uv(32:255))+31
            end
            unique_distb('Tp');                         %Tp
            disp('Done Tp')
            if min(Tp_uv(1:64))~=0
                disp('Unique Tp met')
            else Tp_not_met=find(~Tp_uv(1:64))
            end
            unique_distb('Gp');                         %Gp   
            disp('Done Gp')
            if min(Gp_uv(1:32))~=0
                disp('Unique Gp met')
            elseif looplimits(20)==1 %vsdl
                Gp_not_met=find(~Gp_uv(1:32))
            else disp('No Gp for ADSL')
            end
            unique_distb('Mp');                         %Mp
            disp('Done Mp')
            if min(Mp_uv(1:5))~=0
                disp('Unique Mp met')
            else Mp_not_met=Mp_all(find(~Mp_uv(1:5)))
            end
            Redundancy(looplimits(20));                 %check for redundancy, optimize, print in excel.
        otherwise
    end

else                            % not unique cases
    checkboxStatus = get(handles.framing_check,'Value');%check if to frame
    if (checkboxStatus)
        looplimits(23)=1;       %framing
        if (looplimits(20)==1)  %vdsl - display output format
            disp('Output format: INPmin-DELmax-L-R-q-D-N-1/s-actINP-actDEL-TDR-ADR-Mp-Tp-gp-msgp-perb-per-OR-SEQ-Opi')
        else disp('Output format: INPmin-DELmax-L-R-q-D-N-1/s-actINP-actDEL-TDR-ADR-Mp-Tp-msgp-per-OR-SEQ')
        end
    else looplimits(23)=0;      %no framing
        if (looplimits(20)==1)  %vdsl - display output format
            disp('Output format: INPmin-DELmax-L-R-q-D-N-1/s-actINP-actDEL-TDR-ADR')
        else disp('Output format: INPmin-DELmax-L-R-q-D-N-1/s-actINP-actDEL-TDR-ADR')
        end
    end

     INPrate0_gui(file_name);   %proceed to INP rate
     disp('done');
end

set(handles.start_button,'Enable','On');    %reactivate the start button
set(handles.dummy, 'String','Done');        %update status

function saveloc_button_Callback(hObject, eventdata, handles)
global file_dir file_name

file_dir=uigetdir;              %visual selection of folder
b='\config_file.txt';           %add the default file name
file_name=strcat(file_dir,b);   %global shout
set(handles.saveloc_edit,'String',file_name)    %make changes in path text box

function force_check_Callback(hObject, eventdata, handles)  %fast path case selection
global looplimits

checkboxStatus = get(handles.force_check,'Value');
if(checkboxStatus)                              %fast path case
    looplimits(21)=1;                           %update looplimits            
    
    set(handles.qmin_static,'Enable','Off');    %deactivate the irrelavant
    set(handles.qmax_static,'Enable','Off');
    set(handles.Qmax_edit,'Enable','Off');
    set(handles.Qmin_edit,'Enable','Off');
    set(handles.rmin_static,'Enable','Off');
    set(handles.rmax_static,'Enable','Off');
    set(handles.Rmax_edit,'Enable','Off');
    set(handles.Rmin_edit,'Enable','Off');

else
    looplimits(21)=2;                           %full path
    
    set(handles.qmin_static,'Enable','On');     %if box is unchecked, text is set to normal
    set(handles.qmax_static,'Enable','On');
    set(handles.rmin_static,'Enable','On');
    set(handles.rmax_static,'Enable','On');
    set(handles.Qmax_edit,'Enable','On');
    set(handles.Qmin_edit,'Enable','On');
    set(handles.Rmax_edit,'Enable','On');
    set(handles.Rmin_edit,'Enable','On');
end

function uniquechkbox_Callback(hObject, eventdata, handles)
%global looplimits 

checkboxStatus = get(handles.uniquechkbox,'Value');

if (checkboxStatus)                             %wanted unique values
    %looplimits(25)=1;                          %corresponding looplimit
    set(handles.unique_popmenu,'Enable','On')   %Inactivate irrelavant
    set(handles.delaymin_edit,'Enable','Off')
    set(handles.delaymax_edit,'Enable','Off')
    set(handles.delaystep_edit,'Enable','Off')
    set(handles.INPmin_edit,'Enable','Off')
    set(handles.INPmax_edit,'Enable','Off')
    set(handles.INPstep_edit,'Enable','Off')
    set(handles.text7,'Enable','Off')
    set(handles.text8,'Enable','Off')
    set(handles.text9,'Enable','Off')
    set(handles.text4,'Enable','Off')
    set(handles.text5,'Enable','Off')
    set(handles.text6,'Enable','Off')
    set(handles.framing_check,'Enable','Off')
    set(handles.force_check,'Enable','Off')
    
else %looplimits(25)=0;                         %normal 
    set(handles.unique_popmenu,'Enable','Off')  %activate the necessary
    set(handles.delaymin_edit,'Enable','On')
    set(handles.delaymax_edit,'Enable','On')
    set(handles.delaystep_edit,'Enable','On')
    set(handles.INPmin_edit,'Enable','On')
    set(handles.INPmax_edit,'Enable','On')
    set(handles.INPstep_edit,'Enable','On')
    set(handles.text7,'Enable','On')
    set(handles.text8,'Enable','On')
    set(handles.text9,'Enable','On')
    set(handles.text4,'Enable','On')
    set(handles.text5,'Enable','On')
    set(handles.text6,'Enable','On')
    set(handles.framing_check,'Enable','On')
    set(handles.force_check,'Enable','On')
    set(handles.unique_popmenu,'Value',1)
end

function Lpmin_edit_Callback(hObject, eventdata, handles)
global looplimits

input = str2num(get(hObject,'String'));             %string input from the Lmin box
input2= str2num(get(handles.Lpmax_edit,'String'));  %string input from the Lmax box

vdsel=looplimits(20);                   %changing the Lp limits corresponding to vdsl plan
if (vdsel==1)
    if (isempty(input))||input<8       %cannot be null input
        set(hObject,'String','8');
    elseif input>30000
        set(hObject,'String','30000');
    elseif input> input2                %Lmin cannot be greater than Lmax
        set(hObject,'String',num2str(input2));
    end
else                                    %changing the Lp limits corresponding to adsl plan
    if (isempty(input))||input<8       %cannot be null input         
        set(hObject,'String','8');
    elseif input>7000
        set(hObject,'String','7000');
    elseif input> input2
        set(hObject,'String',num2str(input2));  %Lmin cannot be greater than Lmax
    end
end

guidata(hObject, handles);

function Lpmax_edit_Callback(hObject, eventdata, handles)
global looplimits

input = str2num(get(hObject,'String'));             %string input from the Lmax box
input2= str2num(get(handles.Lpmin_edit,'String'));  %string input from the Lmin box

vdsel=looplimits(20);
if (vdsel==1)                               %vdsl
    if (isempty(input))||input<8
        set(hObject,'String','8');
    elseif input>30000
        set(hObject,'String','30000');
    elseif input< input2
        set(hObject,'String',num2str(input2));  %Lmax cannot be less than Lmin
    end
else                                        %adsl
    if (isempty(input))||input<8
        set(hObject,'String','8');
    elseif input>7000
        set(hObject,'String','7000');
    elseif input< input2
        set(hObject,'String',num2str(input2));  %Lmax cannot be less than Lmin
    end
end
guidata(hObject, handles);

function INPmin_edit_Callback(hObject, eventdata, handles)

input = str2num(get(hObject,'String'));             %start value of INP
input2= str2num(get(handles.INPmax_edit,'String')); %stop value of INP

if (isempty(input))||input<=0           %ensures a non negative non empty value
    set(hObject,'String','1');          %minimum value
elseif input>16
    set(hObject,'String','16');         %maximum value possible
elseif input> input2
    set(hObject,'String',num2str(input2));   %start value cannot be greater than stop
end
guidata(hObject, handles);

function INPmax_edit_Callback(hObject, eventdata, handles)

input = str2num(get(hObject,'String'));             %stop value of INP
input2= str2num(get(handles.INPmin_edit,'String')); %start value of INP

if (isempty(input))||input<=0       %ensures a non negative non empty value
    set(hObject,'String','1');      %minimum INP possible
elseif input>16
    set(hObject,'String','16');     %maximum INP possible
elseif input< input2
    set(hObject,'String',num2str(input2));   %stop value cannot be less than start
end
guidata(hObject, handles);

function delaymin_edit_Callback(hObject, eventdata, handles)

input = str2num(get(hObject,'String'));                 %delay start value
input2= str2num(get(handles.delaymax_edit,'String'));   %delay stop value


if (isempty(input))||input<=0   %ensures a non negative non empty value
    set(hObject,'String','1');  %minimum possible value
elseif input>64
    set(hObject,'String','64'); %maximum possible value
elseif input> input2
    set(hObject,'String',num2str(input2)); %delay start cannot be greater than stop
end
guidata(hObject, handles);

function delaymax_edit_Callback(hObject, eventdata, handles)

input = str2num(get(hObject,'String'));                 %dealy stop value
input2= str2num(get(handles.delaymin_edit,'String'));   %delay start value

if (isempty(input))||input<=0       %ensures a non negative non empty value
    set(hObject,'String','1');      %minimum possible value
elseif input>64
    set(hObject,'String','64');     %maximum possible value
elseif input< input2
    set(hObject,'String',num2str(input2));   %stop value cannot be less than start
end
guidata(hObject, handles);

function Rmin_edit_Callback(hObject, eventdata, handles)

input = round(str2num(get(hObject,'String')));              %int min range of R
input2= round(str2num(get(handles.Rmax_edit,'String')));    %int max range of R

if mod(input,2)~=0                           %necessarily multiple of 2
    set(hObject,'String',num2str(input+1));  %ensuring multiple of 2
end

if (isempty(input))||input<=0        %non empty non zero value
    set(hObject,'String','2');       %minimum possible value
elseif input>16
    set(hObject,'String','16');      %maximum possible value
elseif input> input2
    set(hObject,'String',num2str(input2)); %ensure min not greater than max
else set(hObject,'String',num2str(input));
end
guidata(hObject, handles);

function Rmax_edit_Callback(hObject, eventdata, handles)

input = round(str2num(get(hObject,'String')));              %int max range of R
input2= round(str2num(get(handles.Rmin_edit,'String')));    %int min range of R

if mod(input,2)~=0                          %necessarily multiple of 2
    set(hObject,'String',num2str(input+1));
end

if (isempty(input))||input<=0               %non empty non zero value
    set(hObject,'String','2');              %minimum possible
elseif input>16
    set(hObject,'String','16');             %maximum possible
elseif input< input2
    set(hObject,'String',num2str(input2));  %max cannot be less than mim
else set(hObject,'String',num2str(input));
end
guidata(hObject, handles);

function Qmin_edit_Callback(hObject, eventdata, handles)

input = round(str2num(get(hObject,'String')));              %minimum possible value of q
input2= round(str2num(get(handles.Qmax_edit,'String')));    %maximum possible value of q

if input>input2
    set(hObject,'String',num2str(input2));      %ensure qmax, qmin
elseif input < 1
    set(hObject,'String','1');                  %cannot be less than 1
else set(hObject,'String',num2str(input));
end

function Qmax_edit_Callback(hObject, eventdata, handles)

input = round(str2num(get(hObject,'String')));          %max possible q value
input2= round(str2num(get(handles.Qmin_edit,'String')));%min possible q value
input2 = ceil(input2);
if input<input2
    set(hObject,'String',num2str(input2));  %ensure qmax, qmin
elseif input > 8 
    set(hObject,'String','8');              %ensure qmax not >8
else set(hObject,'String',num2str(input));
end

function msg_edit_Callback(hObject, eventdata, handles)

global looplimits

input = str2num(get(hObject,'String')); %value entered
vdsel=looplimits(20);                   %vdsl selection

if (vdsel==1)                           %vdsl
    if (isempty(input))||input<16
        set(hObject,'String','16');     %min limits
    elseif input>248
        set(hObject,'String','248');    %maximum limit
    end
else                                    %adsl
    if (isempty(input))||input<4
        set(hObject,'String','4');      %min limts 
    elseif input>64
        set(hObject,'String','64');     %max limits      
    end
end

guidata(hObject, handles);


function plan_set_menu_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function msg_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function unique_popmenu_Callback(hObject, eventdata, handles)
function unique_popmenu_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function varargout = gui_inp_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;
function Qmax_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function saveloc_edit_Callback(hObject, eventdata, handles)
function framing_check_Callback(hObject, eventdata, handles)
function saveloc_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Qmin_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Rmax_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Rmin_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function INPmax_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function INPstep_edit_Callback(hObject, eventdata, handles)
function INPstep_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Lpmax_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Lpstep_edit_Callback(hObject, eventdata, handles)
function Lpstep_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function INPmin_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function delaymax_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function delaystep_edit_Callback(hObject, eventdata, handles)
function delaystep_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Lpmin_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function delaymin_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

