function unique_inp1(var_name,inp)
%inputs are var_name which has to be unique value. possible values are
%[D,N,Tp,Gp and Mp] inp is the INP value to be regarded. This function
%necessarily needs unique_frame file to complete.

global looplimits 
global D_uv Gp_uv Tp_uv N_uv Mp_uv inp_co inp_co_uv

INP_start   = looplimits(1);    %user defined
INP_end     = looplimits(2);
INP_step    = looplimits(3);

L_min       = looplimits(7);    %user defined
L_max       = looplimits(8);
L_step      = looplimits(9);

R_min       = looplimits(10);   %user defined
R_max       = looplimits(11);

Q_min       = looplimits(12);   %user defined
Q_max       = looplimits(13);

D_max       = looplimits(14);   %profile defined

inv_S_min   = looplimits(15);   %profile defined
S_min       = 1.0/inv_S_min;
S_max       = looplimits(16);

IntlMem_max = looplimits(17);   %profile defined (+user defined %)

Fs          = looplimits(18);   %profile defined
msg_min     = looplimits(19);   %user defined

if looplimits(20)==1    %adsl/vdsl
    sl='vdsl';
elseif looplimits(20)==2
    sl='adsl';
end

if looplimits(22)==1    %possible values of D in adsl
    dlim=21;            %downstream
else dlim=7;            %upstream
end
n_chk=0;
D_adsl=[1,2,4,8,16,32,64,96,128,160,192,224,256,288,320,352,384,416,448,480,511];   %Possible D values adsl.

IntlMem_max=round(IntlMem_max*.80);

if strcmp('vdsl',sl)
    if strcmp('D',var_name)                 %unique D value
        for Lp=L_max:-L_step:L_min
             Lp1=Lp+randint(1); 
            
            for n=255:-1:32
                if  min(D_uv(1:D_max))==0   %return if all D met
                    s=((8*n)/Lp1);           %Table 9-6/G.993.2 Derived characteristics 
                    if s<=S_max && s>=S_min
                        %r=round((-14*length(find(D_uv(1:D_max)))+49136)/3071);
                        r=round((-14*length(find(D_uv(1:D_max)))+(16*(D_max-1)))/(D_max-1));
                        r=min(16,max(2,r+(mod(r,2)~=0)));
                        for q=Q_min:min(r/2,Q_max)  %avoiding division by zero.
                                d=ceil(inp*n/(floor(r/(2*q))*s));   
                                if d<=D_max
                                    if D_uv(d)==0       %check if this has previously occured; unique D?
                                        if ((n/q)-1)*(d-1)< (IntlMem_max+(length(find(D_uv(1:D_max)))*(1-IntlMem_max/(D_max-1))))     %Table 9-3/G.993.2
                                            if (inp_co_uv(inp)<inp_co(inp)) %proceed iff INP(count<desired_count)
                                                s_temp1 = d;                %co-prime criteria of n and i
                                                s_temp2 = (n/q);            %Table 9-3/G.993.2
                                                while s_temp2 > 0           %unless remainder = 0
                                                    s_temp3 = floor(s_temp1 / s_temp2);
                                                    s_temp4 = s_temp1 - s_temp2 * s_temp3;
                                                    s_temp1 = s_temp2;
                                                    s_temp2 = s_temp4;
                                                end
                                                if s_temp1 == 1             %valid co-primes

                                                    frame(1)=inp;           %proceed to framing
                                                    frame(2)=Lp1;
                                                    frame(3)=n;
                                                    frame(4)=s;
                                                    frame(5)=d;
                                                    frame(6)=q;
                                                    frame(7)=r;
                                                 
                                                    unique_frame1(frame,var_name,sl,msg_min);
                                                end
                                            else return     %INP count met, next INP from unique_distrb
                                            end
                                            
                                        end
                                    end
                                end
                            end
                        %end
                    end
                else return;
                end
            end
        end
    elseif strcmp('N',var_name)                         %unique N value IntlMem_max+(length(find(N_uv(32:255)))*(1-IntlMem_max/223)
        for Lp=L_max:-L_step:L_min
            Lp1=Lp+randint(1);
            for n=255:-1:32
                if  min(N_uv(32:255))==0                %return if all N met
                    s=((8*n)/Lp1);
                    if s<=S_max && s>=S_min
                        r=round((-14*length(find(N_uv(32:255)))+3568)/223);
                        r=min(16,max(2,r+(mod(r,2)~=0)));
                        %for r=R_max:-2:R_min
                            for q=Q_min:min(r/2,Q_max)  %avoiding division by 0
                                d=ceil(inp*n/(floor(r/(2*q))*s));
                                if d<=D_max && d>1
                                    if N_uv(n)==0       %check and skip if it has already occured
                                        if ((n/q)-1)*(d-1)< (IntlMem_max+(length(find(N_uv(32:255)))*(1-IntlMem_max/223)))%(-(length(find(N_uv(32:255)))*IntlMem_max)/223+98300)    %maximum interleaver
                                            if (inp_co_uv(inp)<inp_co(inp)) %proceed iff INP(count<desired_count)
                                                s_temp1 = d;                %co-prime criteria of n and i
                                                s_temp2 = (n/q);
                                                while s_temp2 > 0           %unless remainder = 0
                                                    s_temp3 = floor(s_temp1 / s_temp2);
                                                    s_temp4 = s_temp1 - s_temp2 * s_temp3;
                                                    s_temp1 = s_temp2;
                                                    s_temp2 = s_temp4;
                                                end
                                                if s_temp1 == 1             %valid co-prime

                                                    frame(1)=inp;
                                                    frame(2)=Lp1;
                                                    frame(3)=n;
                                                    frame(4)=s;
                                                    frame(5)=d;
                                                    frame(6)=q;
                                                    frame(7)=r;             %proceed to framing
                                              
                                                    unique_frame1(frame,var_name,sl,msg_min);
                                                end
                                            else return
                                            end
                                        end
                                     end
                                end
                            end
                        %end
                    end
                else return;
                end
            end

        end

    elseif strcmp('Tp',var_name)        %unique cases in Tp
        for Lp=L_max:-L_step:L_min
            for n=255:-1:32
                 Lp1=Lp+randint(1);
                if  min(Tp_uv(1:64))==0 %return if all Tp met
                    s=((8*n)/Lp1);
                    if s<=S_max && s>=S_min
                        %for r=R_max:-2:R_min
                        r=round((-2*length(find(Tp_uv(1:64)))+144)/9);
                        r=min(16,max(2,r+(mod(r,2)~=0)));    
                        for q=Q_min:min(r/2,Q_max)
                                d=ceil(inp*n/(floor(r/(2*q))*s));
                                if d<=D_max && d>1
                                    if ((n/q)-1)*(d-1)< (IntlMem_max+(length(find(Tp_uv(1:64)))*(1-IntlMem_max/63)))     %maximum interleaver
                                        if (inp_co_uv(inp)<inp_co(inp)) %proceed iff INP(count<desired_count)
                                            s_temp1 = d;                %co-prime criteria of n and i
                                            s_temp2 = (n/q);
                                            while s_temp2 > 0           %unless remainder = 0
                                                s_temp3 = floor(s_temp1 / s_temp2);
                                                s_temp4 = s_temp1 - s_temp2 * s_temp3;
                                                s_temp1 = s_temp2;
                                                s_temp2 = s_temp4;
                                            end
                                            if s_temp1 == 1             %valid co-prime

                                                frame(1)=inp;
                                                frame(2)=Lp1;
                                                frame(3)=n;
                                                frame(4)=s;
                                                frame(5)=d;
                                                frame(6)=q;
                                                frame(7)=r;
                                                unique_frame1(frame,var_name,sl,msg_min);
                                            end
                                        else return
                                        end
     
                                    end
                                end
                            end
                        %end
                    end
                else return;
                end
            end
        end

    elseif strcmp('Gp',var_name)
        for Lp=L_min:(L_step):L_max
            for n=255:-1:32
                if  min(Gp_uv(1:32))==0 %return if all possible Gp met
                    s=((8*n)/Lp);
                    if s<=S_max && s>=S_min
                        for r=R_max:-2:R_min
                            for q=Q_min:min(r/2,Q_max)
                                d=ceil(inp*n/(floor(r/(2*q))*s));
                                if d<=D_max && d>1
                                    if ((n/q)-1)*(d-1)< (IntlMem_max+(length(find(Gp_uv(1:32)))*(1-IntlMem_max/31)))     %maximum interleaver case
                                        if (inp_co_uv(inp)<inp_co(inp)) %proceed iff INP(count<desired_count)
                                            s_temp1 = d;                %co-prime criteria of n and i
                                            s_temp2 = (n/q);
                                            while s_temp2 > 0           %unless remainder = 0
                                                s_temp3 = floor(s_temp1 / s_temp2);
                                                s_temp4 = s_temp1 - s_temp2 * s_temp3;
                                                s_temp1 = s_temp2;
                                                s_temp2 = s_temp4;
                                            end
                                            if s_temp1 == 1             %valid co-prime

                                                frame(1)=inp;
                                                frame(2)=Lp;
                                                frame(3)=n;
                                                frame(4)=s;
                                                frame(5)=d;
                                                frame(6)=q;
                                                frame(7)=r;             %proceed to framing
                                                unique_frame1(frame,var_name,sl,msg_min);
                                            end
                                        else return
                                        end
                                    end
                                end
                            end
                           
                        end
                    end
                else return;
                end
            end
        end
    elseif strcmp('Mp',var_name)
        for Lp=L_min:L_step:L_max
            for n=255:-1:32
                if  min(Mp_uv(1:5))==0  %check and return if all possible Mp are met.
                    s=((8*n)/Lp);
                    if s<=S_max && s>=S_min
                        for r=R_max:-2:R_min
                            for q=Q_min:min(r/2,Q_max)
                                d=ceil(inp*n/(floor(r/(2*q))*s));
                                if d<=D_max && d>1
                                    if ((n/q)-1)*(d-1)< (IntlMem_max+(length(find(Mp_uv(1:5)))*(1-IntlMem_max/4)))
                                        if (inp_co_uv(inp)<inp_co(inp)) %proceed iff INP(count<desired_count)
                                            s_temp1 = d;                %co-prime criteria of n and i
                                            s_temp2 = (n/q);
                                            while s_temp2 > 0           %unless remainder = 0
                                                s_temp3 = floor(s_temp1 / s_temp2);
                                                s_temp4 = s_temp1 - s_temp2 * s_temp3;
                                                s_temp1 = s_temp2;
                                                s_temp2 = s_temp4;
                                            end
                                            if s_temp1 == 1             %valid co-prime

                                                frame(1)=inp;
                                                frame(2)=Lp;
                                                frame(3)=n;
                                                frame(4)=s;
                                                frame(5)=d;
                                                frame(6)=q;
                                                frame(7)=r;             %proceed to framing
                                                unique_frame1(frame,var_name,sl,msg_min);
                                            end
                                        else return
                                        end
                                    end
                                end
                            end
                        end
                       
                    end
                else return;
                end
            end
        end

    end
elseif strcmp('adsl',sl)                        %adsl cases
    if strcmp('D',var_name)
        for Lp=L_max:-L_step:L_min
            for n=255:-1:32
%                 if mod(n,2)==0 
%                     n1=n+1;
%                     n_chk=1;
%                 else n_chk=0;
%                     n1=n;
%                 end
                if  min(D_uv(1:dlim))==0        %all possible D values met?
                    s=((8*n)/Lp);               %Table 7-7/G.992.3
                    if s<=S_max && s>=S_min
                        for r=R_max:-2:R_min
                            q=1;                %adsl no q(=1)
                            d=ceil(inp*n/(floor(r/(2*q))*s));
                            d=getD_adsl(d,looplimits(22));              %ceil to the next nearest valid D
                            if d<=D_max
                                if (D_uv(find(D_adsl==d))==0)           %iff this D has not occured previously
                                    if ((n/q)-1)*(d-1)< IntlMem_max     %maximum interleaver case
                                        if (inp_co_uv(inp)<inp_co(inp)) %proceed iff INP(count<desired_count)
                                            s_temp1 = d;                %co-prime criteria of n and i
                                            s_temp2 = (n/q);
                                            while s_temp2 > 0           %unless remainder = 0
                                                s_temp3 = floor(s_temp1 / s_temp2);
                                                s_temp4 = s_temp1 - s_temp2 * s_temp3;
                                                s_temp1 = s_temp2;
                                                s_temp2 = s_temp4;
                                            end
                                            if s_temp1 == 1             %valid co-prime

                                                frame(1)=inp;
                                                frame(2)=Lp;
                                                frame(3)=n;%n1-n_chk;
                                                frame(4)=s;
                                                frame(5)=d;
                                                frame(6)=q;
                                                frame(7)=r;             %proceed to framing
                                                
                                                unique_frame1(frame,var_name,sl,msg_min);

                                            end
                                        else return
                                        end
                                    end
                                end
                            end
                          
                        end
                    end
                else return;
                end
            end
        end
    elseif strcmp('N',var_name)
        for Lp=L_max:-L_step:L_min
            for n=255:-1:32
                if  min(N_uv(32:255))==0
                    if  N_uv(n)==0
                        if mod(n,2)==0
                            n1=n+1;
                            n_chk=1;
                        else n_chk=0;
                            n1=n;
                        end
                        s=((8*n1)/Lp);
                        if s<=S_max && s>=S_min
                            r=round((-14*length(find(N_uv(32:255)))+3568)/223);
                        r=min(16,max(2,r+(mod(r,2)~=0)));
                            %for r=R_max:-2:R_min
                                q=1;
                                d=ceil(inp*n1/(floor(r/(2*q))*s));
                                d=getD_adsl(d,looplimits(22));
                                if d<=D_max

                                    if ((n1/q)-1)*(d-1)< (IntlMem_max+(length(find(N_uv(32:255)))*(1-IntlMem_max/223)))
                                        if (inp_co_uv(inp)<inp_co(inp))
                                            s_temp1 = d;            %co-prime criteria of n and i
                                            s_temp2 = (n1/q);        %proceed iff INP(count<desired_count)
                                            while s_temp2 > 0       %unless remainder = 0
                                                s_temp3 = floor(s_temp1 / s_temp2);
                                                s_temp4 = s_temp1 - s_temp2 * s_temp3;
                                                s_temp1 = s_temp2;
                                                s_temp2 = s_temp4;
                                            end
                                            if s_temp1 == 1         %valid co-prime

                                                frame(1)=inp;
                                                frame(2)=Lp;
                                                frame(3)=n1-n_chk;
                                                frame(4)=s;
                                                frame(5)=d;
                                                frame(6)=q;
                                                frame(7)=r;         %proceed to framing

                                                unique_frame1(frame,var_name,sl,msg_min);
                                            end
                                        else return
                                        end
                                    end
                                end
                            %end
                        end
                    end
                else return;
                end
            end
        end

    elseif strcmp('Tp',var_name)
        for Lp=L_max:-L_step:L_min
            for n=255:-1:32
%                 if mod(n,2)==0
%                     n1=n+1;
%                     n_chk=1;
%                 else n_chk=0;
%                     n1=n;
%                 end
                if  min(Tp_uv(1:64))==0     %check if all Tp have occured
                    s=((8*n)/Lp);
                    if s<=S_max && s>=S_min
                        r=round((-2*length(find(Tp_uv(1:64)))+144)/9);
                        r=min(16,max(2,r+(mod(r,2)~=0)));
                        %for r=R_max:-2:R_min
                            q=1;
                            d=ceil(inp*n/(floor(r/(2*q))*s));
                            d=getD_adsl(d,looplimits(22));
                            if d<=D_max
                                if ((n/q)-1)*(d-1)< (IntlMem_max+(length(find(Tp_uv(1:64)))*(1-IntlMem_max/63)) )
                                    if (inp_co_uv(inp)<inp_co(inp))     %proceed iff INP(count<desired_count)
                                        s_temp1 = d;                    %co-prime criteria of n and i
                                        s_temp2 = (n/q);

                                        while s_temp2 > 0               %unless remainder = 0
                                            s_temp3 = floor(s_temp1 / s_temp2);
                                            s_temp4 = s_temp1 - s_temp2 * s_temp3;
                                            s_temp1 = s_temp2;
                                            s_temp2 = s_temp4;
                                        end

                                        if s_temp1 == 1         %valid co-prime

                                            frame(1)=inp;
                                            frame(2)=Lp;
                                            frame(3)=n;%1-n_chk;
                                            frame(4)=s;
                                            frame(5)=d;
                                            frame(6)=q;
                                            frame(7)=r;         %proceed to framing
                                            
                                            unique_frame1(frame,var_name,sl,msg_min);
                                        end
                                    else return
                                    end
                                end
                            end
                        %end
                    end
                else return;
                end
            end
        end

    elseif strcmp('Mp',var_name)
        for Lp=L_max:-L_step:L_min
            for n=255:-1:32
              
                if  min(Mp_uv(1:5))==0      %check that all Mp have occured
                    s=((8*n)/Lp);
                    if s<=S_max && s>=S_min
                        for r=R_max:-2:R_min
                            q=1;
                            d=ceil(inp*n/(floor(r/(2*q))*s));
                            d=getD_adsl(d,looplimits(22));
                            if d<=D_max
                                if ((n/q)-1)*(d-1)< (IntlMem_max+(length(find(Mp_uv(1:5)))*(1-IntlMem_max/4)))
                                    if (inp_co_uv(inp)<inp_co(inp))     %proceed iff INP(count<desired_count)
                                        s_temp1 = d;                    %co-prime criteria of n and i
                                        s_temp2 = (n/q);
                                        while s_temp2 > 0               %unless remainder = 0
                                            s_temp3 = floor(s_temp1 / s_temp2);
                                            s_temp4 = s_temp1 - s_temp2 * s_temp3;
                                            s_temp1 = s_temp2;
                                            s_temp2 = s_temp4;
                                        end
                                        if s_temp1 == 1                 %valid co-prime

                                            frame(1)=inp;
                                            frame(2)=Lp;
                                            frame(3)=n;%n1-n_chk;
                                            frame(4)=s;
                                            frame(5)=d;
                                            frame(6)=q;
                                            frame(7)=r;                 %proceed to framing
                                            
                                            unique_frame1(frame,var_name,sl,msg_min);
                                        end
                                    else return
                                    end
                                end
                            end
                        end
                        
                    end
                else return;
                end
            end
        end
    end
end

function [d]=getD_adsl(d,stream)
if stream==1 %downstream, Table 7-8/G.992.3
    if d==1 ||d==2
        d;
    elseif d>2 && d<=4
        d=4;
    elseif d>4 && d<=8
        d=8;
    elseif d>8 && d<=16
        d=16;
    elseif d>16 && d<=32
        d=32;
    elseif d>32 && d<=64
        d=64;
    elseif d>64 && d<=96
        d=96;
    elseif d>96 && d<=128
        d=128;
    elseif d>128 && d<=160
        d=160;
    elseif d>160 && d<=192
        d=192;
    elseif d>192 && d<=224
        d=224;
    elseif d>224 && d<=256
        d=256;
    elseif d>256 && d<=288
        d=288;
    elseif d>288 && d<=320
        d=320;
    elseif d>320 && d<=352
        d=352;
    elseif d>352 && d<=384
        d=384;
    elseif d>384 && d<=416
        d=416;
    elseif d>416 && d<=448
        d=448;
    elseif d>448 && d<=480
        d=480;
    elseif d>480 && d<=511
        d=511;
    end
elseif stream==2 %upstream
    if d==1 || d==2
        d;
    elseif d>2 && d<=4
        d=4;
    elseif d>4 && d<=8
        d=8;
    elseif d>8 && d<=16
        d=16;
    elseif d>16 && d<=32
        d=32;
    elseif d>32 && d<=64
        d=64;
    end
end



