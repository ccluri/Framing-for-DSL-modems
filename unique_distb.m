function unique_distb(var_name)
%var_name represents the parameter intended to find the unique values.
%possible values var_name takes (D,N,Tp,Gp,Mp), This funtion calls for
%another funtion unique_INP

global looplimits
global inp_co inp_co_uv inp_co_def

%inp_co_def - proportions as much defined 
%inp_co - as much it has to be
%inp_co_uv - as much as it

mx_inp=16; %maximum inp set to 16

if strcmp('D',var_name)
    inp_co_uv(1:mx_inp)=0;
    if looplimits(20)==1        %vdsl
        inp_co(1:mx_inp)=round(inp_co_def(1:mx_inp).*looplimits(14));
    elseif looplimits(22)==1    %adsl ds
        inp_co(1:mx_inp)=round(inp_co_def(1:mx_inp).*21);
    elseif looplimits(22)==2    %adsl us
        inp_co(1:mx_inp)=round(inp_co_def(1:mx_inp).*7);
    end
elseif strcmp('N',var_name)
    inp_co_uv(1:mx_inp)=0;
    inp_co(1:mx_inp)=round(inp_co_def(1:mx_inp).*255);  %32-255
elseif strcmp('Tp',var_name)
    inp_co_uv(1:mx_inp)=0;
    inp_co(1:mx_inp)=round(inp_co_def(1:mx_inp).*64);   %1-64
elseif strcmp('Gp',var_name)
    inp_co_uv(1:mx_inp)=0;
    inp_co(1:mx_inp)=round(inp_co_def(1:mx_inp).*32);   %1-32
elseif strcmp('Mp',var_name)
    inp_co_uv(1:mx_inp)=0;
    inp_co(1:mx_inp)=round(inp_co_def(1:mx_inp).*5);    %[5]
end

for inp=2:mx_inp;               %maximum number of cases with INP=2:16 (waterfilling algorithm)
    unique_inp1(var_name,inp);  %distributed function
end

disp('Level-1-1');                %status print

inp_co(2:16)=inp_co(2:16).*8;   %resetting maximum # of INP cases, inp=2:16 leak proof check & waterfilling INP=1.
for inp=1:16
unique_inp1(var_name,inp);      %distributed function    
end

disp('Level-1-2');                %status print

inp_co(1)=inp_co(1).*8;         %leak proof check INP=1
unique_inp1(var_name,1);        %distributed



inp_co(1:16)=inp_co(1:16)./8;   %RESETTING, distrb.

for inp=2:mx_inp;               %maximum number of cases with INP=2:16 (waterfilling algorithm)
    unique_inp(var_name,inp);
end

disp('Level-2-1');                %status print

inp_co(2:16)=inp_co(2:16).*8;   %resetting maximum # of INP cases, inp=2:16 leak proof check & waterfilling INP=1.
for inp=1:16
unique_inp(var_name,inp);
end

disp('Level-2-2');                %status print

inp_co(1)=inp_co(1).*8;         %leak proof check INP=1
unique_inp(var_name,1);

return;

