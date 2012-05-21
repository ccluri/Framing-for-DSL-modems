function frame_gui(Fs,sl,msgp_min,output,fid)
% this program, checks if a particular case is framable based on output
% configarations form INPrate_gui.m and if it is framable, it
% saves the framing parameters into the file_name inputed from the
% INPrate_gui.m file
%
% the cases of adsl and vdsl are treated separately, output variable is the
% input, sl represents the choice of adsl or vdsl, msg_min is userdefined
% input, fid corresponds is the handle to the file_name.txt

% Algorithm for framing -vdsl
% the Q1 for perb calculations based on Tdr is evaluated
% the validity of framin is checked while running over limits of Tp(1-64),
% over all possible Mp values.
% per is between [1&20] and Tp/Mp is an integer
% Gp lies between [1&32]
% Mp/s<64 rule1 and
% (floor(Gp/Tp)*d)+(floor(d/Tp)*(c))+min(mod(d,Tp),c)<=8 rule2

% Algorithm for framing -adsl
% the validity of framin is checked while running over limits of Tp(1-64),
% over all possible Mp values.
% OR=(Mp*32*a)/(Tp)
% OR lies between[0.8&64]
% and per lies between [15&20]

Ls  =   output(3);
a   =   output(8);
r   =   output(4);
nfec=   output(7);
s   =   1/a;

Mp_all=[1,2,4,8,16];

chk=0;

if strcmp('vdsl',sl)    %check if its vdsl

    op(1:32)=0;     %default values
    %perb calculations -Table 9-6/G.993.2,,Derived framing parameters 
    tdr=Ls*Fs*10^-3; %kbits/sec
    if tdr>=7880
        q1=17000;
    else q1=17000*(tdr/7880);
    end

    for Tp=1:64
        for i=1:5
            Mp=Mp_all(i);
            perb=((Tp*nfec)/Mp)*floor((q1*Mp)/(Tp*nfec));   %bytes
            per= (8*perb*1000)/(Ls*Fs);                     %Table 9-6/G.993.2,Derived framing parameters
            if per<=20 && per>=1 && mod(Tp,Mp)==0           %per lies between 1ms and 20ms
                msgp=msgp_min;
                Gp=ceil(((msgp*Tp*1000)/(8*a*Mp*Fs))+((6*nfec*Tp)/(perb*Mp)));  %OR and msgp relations(Table 9-6/G.993.2) 
                if Gp<=32 && Gp>=1      %Table 9-6/G.993.2,Primary framing parameters
                    
                    if (Mp*a)<=64               %rule1, 9.5.2.1/G.993.2
                        d=(ceil(Mp*a));
                        c=mod(Gp,Tp);

                        r2v1=(floor(Gp/Tp)*d);  %rule2, 9.5.2.1/G.993.2
                        r2v2= (floor(d/Tp)*(c));
                        r2v3=min(mod(d,Tp),c);
                        r2vf=r2v1+r2v2+r2v3;

                        if r2vf <=8                      %rule2 (floor(Gp/Tp)*d)+(floor(d/Tp)*(c))+min(mod(d,Tp),c)
                            if mod((nfec-r),Mp)==0       %ensuring that Bpn-Nfec integer conditions are met. (Figure 9-2/G.993.2)
                                chk=1;                   %check flag
                                op(1:32)=0;              %initialize/clear
                                for x=1:Tp
                                    if x<=((Gp)-(Tp*floor(Gp/Tp)))  %9.5.2.1/G.993.2
                                        op(x)=ceil(Gp/Tp);
                                    else op(x)=floor(Gp/Tp);
                                    end
                                end
                            end
                        end

                    end

                    

                    if chk==1
                        chk=0;

                        OR=((Gp*Mp*8*Fs*a)/(Tp*1000)); %Table 9-6/G.993.2,Derived framing parameters
                        SEQ=((perb*Mp*Gp)/(nfec*Tp));
                        msgp=OR*(((SEQ)-6)/(SEQ));

                        output(13)=Mp;
                        output(14)=Tp;
                        output(15)=Gp;
                        output(16)=msgp;
                        output(17)=perb;
                        output(18)=per;
                        output(19)=OR;
                        output(20)=SEQ;

                        output(21:52)=op(1:32);
                        fprintf(fid,'%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g',output(1:52));
                        fprintf(fid,'\n');
                    end
                end
            end
        end
    end
end


if strcmp('adsl',sl)
    for Tp=1:64
        for i=1:5
            Mp=Mp_all(i);
            if mod((nfec-r-ceil(Mp/Tp)),Mp)==0  %ensuring Bpn-nfec interger condition, Figure 7-7/G.992.3
                if (s>=(Mp/16))&&(s<=(32*Mp))   %Sp-Mp conditon, Table 7-8/G.992.3
                    if r==0                     %fast path case
                        Mp=1;                   %Table 7-8/G.992.3 – Valid framing configurations
                    end

                    OR=(Mp*32*a)/(Tp);                      %kbit/s, Table 7-7/G.992.3
                    if OR>=0.8 && OR<=64 && OR >msgp_min    %Overhead Rate Constraints, Table 7-8/G.992.3
                        msgc=ceil((6*msgp_min)/((OR)-msgp_min));    %msgp=MSGmin, using rel from Framing_VDSL_V05f(version1).xls adsl2
                        SEQ=msgc+6;                         %Table 7-7/G.992.3
                        per= (Tp*SEQ)/(4*a*Mp);             %Table 7-7/G.992.3

                        if per<=20 && per>=15           %Table 7-8/G.992.3,Overhead Channel Period

                            output(13)=Mp;
                            output(14)=Tp;
                            output(15)=msgc;
                            output(16)=per;
                            output(17)=OR;
                            output(18)=SEQ;

                            fprintf(fid,'%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g',output(1:18));
                            fprintf(fid,'\n');
                        end
                    end
                end
            end
        end
    end
end


