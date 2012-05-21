function unique_frame1(frame,var_name,sl,msgp)
%inputs for unique_frame are frame(from unique_gui); var_name is the 
%parameter to which unique values aer to be found (D,N,Tp,Gp,Mp); sl, 
%is the string that dictates if vdsl/adsl, msgp is the MSGmin rate.
%the value returned is a text file, the values correspond to that as in
%output line printed
global D_uv Gp_uv Tp_uv N_uv Mp_uv inp_co inp_co_uv file_dir
chk=0;

Lp=frame(2);    %from unique_gui
n=frame(3);
s=frame(4);
a=1/s;
d=frame(5);
r=frame(7);

D_max=2048; % ****************************


Mp_all=[1,2,4,8,16];        %all possible Mp values, Table 7-8/G.992.3 & Table 9-6/G.993.2
D_adsl=[1,2,4,8,16,32,64,96,128,160,192,224,256,288,320,352,384,416,448,480,511];   %adsl possible D values,Table 7-8/G.992.3

tdr=Lp*4000*10^-3;          %kbits/sec; Lp=((8*n)/s);
if tdr>=7880                %Table 9-6/G.993.2,Derived framing parameters
    q1=17000;
else q1=17000*(tdr/7880);
end
%op(1:32)=[];
if strcmp('vdsl',sl)
    if strcmp('D',var_name)     %unique values of D,vdsl
        for Tp=1:64             %possible Tp vals
            %for i=1:5           %pointer to Mp vals
             i=max(1,min(5,round((4*length(find(D_uv(1:D_max)))+(D_max-1))/D_max)));
                Mp=Mp_all(i);
                perb=((Tp*n)/Mp)*floor((q1*Mp)/(Tp*n)); %bytes
                per= (8*perb*1000)/(Lp*4000);           %Table 9-6/G.993.2,Derived framing parameters
                if per<=20 && per>=1 && mod(Tp,Mp)==0   %per lies between 1ms and 20ms
                    Gp=ceil(((msgp*Tp*s*1000)/(8*Mp*4000))+((6*n*Tp)/(perb*Mp)));   %OR and msgp relations(Table 9-6/G.993.2) 
                    if Gp<=32 && Gp>=1                  %Table 9-6/G.993.2,Primary framing parameters
                        if (Mp/s)<=64                   %rule1, 9.5.2.1/G.993.2
                            de=(ceil(Mp/s));
                            ce=mod(Gp,Tp);              %rule2, 9.5.2.1/G.993.2
                            if (floor(Gp/Tp)*de)+(floor(de/Tp)*(ce))+min(mod(de,Tp),ce)<=8  
                                if D_uv(d)==0           %Is this D value unique?
                                    if mod((n-r),Mp)==0 %ensuring that Bpn-Nfec integer conditions are met. (Figure 9-2/G.993.2)
                                        chk=1;          %flag for update
                                        op(1:32)=0;     %initialize/clear 
                                        for x=1:Tp
                                            if x<=((Gp)-(Tp*floor(Gp/Tp)))  %9.5.2.1/G.993.2
                                                op(x)=ceil(Gp/Tp);
                                            else op(x)=floor(Gp/Tp);
                                            end
                                        end
                                    end
                                else return;
                                end
                            end
                        end

                        if chk==1
                            chk=0;

                            Gp_uv(Gp)=1;    %flag Gp
                            Tp_uv(Tp)=1;    %flag Tp
                            D_uv(d)=1;      %flag D
                            N_uv(n)=1;      %flag N
                            Mp_uv(i)=1;     %flag Mp
                            inp_co_uv(frame(1))=inp_co_uv(frame(1))+1; %incrementing the occurance of an INP

                            frame(8)=msgp;
                            frame(9)=Tp;
                            frame(10)=Mp;
                            frame(11)=Gp;
                            frame(12)=perb;
                            frame(13)=per;
                            frame(14)=(perb*Mp*Gp)/(n*Tp);  %SEQ
                            frame(15:46)=op(1:32);

                            dlmwrite(strcat(file_dir,'\Unique-D-vdsl.txt'), frame, 'newline','pc','-append')
                        end
                    end
                end
           % end
        end
    elseif strcmp('N',var_name) %unique in N values
        %Tp=min(64,max(1,round((-63*length(find(N_uv(32:255)))+14335)/223)));
        for Tp=1:64             %possible Tp vals
            %for i=1:5           %pointer to Mp vals
            i=max(1,min(5,round((4*length(find(N_uv(32:255)))+223)/223)));
                Mp=Mp_all(i);
                perb=((Tp*n)/Mp)*floor((q1*Mp)/(Tp*n)); %bytes
                per= (8*perb*1000)/(Lp*4000);
                if per<=20 && per>=1 && mod(Tp,Mp)==0   %per lies between 1ms and 20ms
                    Gp=ceil(((msgp*Tp*s*1000)/(8*Mp*4000))+((6*n*Tp)/(perb*Mp)));   %msgp,OR relation
                    if Gp<=32 && Gp>=1
                        if (Mp/s)<=64                   %rule1
                            de=(ceil(Mp/s));
                            ce=mod(Gp,Tp);              %rule2
                            if (floor(Gp/Tp)*de)+(floor(de/Tp)*(ce))+min(mod(de,Tp),ce)<=8 
                                if mod((n-r),Mp)==0     %bpn-nfec integer
                                    if N_uv(n)==0       %is this a unique occurance of N?
                                        chk=1;
                                        op(1:32)=0;
                                        for x=1:Tp
                                            if x<=((Gp)-(Tp*floor(Gp/Tp)))
                                                op(x)=ceil(Gp/Tp);
                                            else op(x)=floor(Gp/Tp);
                                            end
                                        end
                                    else return;
                                    end
                                end
                            end
                        end

                        if chk==1
                            chk=0;
                            Gp_uv(Gp)=1;    %flags
                            Tp_uv(Tp)=1;
                            D_uv(d)=1;
                            N_uv(n)=1;
                            Mp_uv(i)=1;
                            inp_co_uv(frame(1))=inp_co_uv(frame(1))+1; %increase INP counts

                            frame(8)=msgp;
                            frame(9)=Tp;
                            frame(10)=Mp;
                            frame(11)=Gp;
                            frame(12)=perb;
                            frame(13)=per;
                            frame(14)=(perb*Mp*Gp)/(n*Tp);  %SEQ
                            frame(15:46)=op(1:32);

                            dlmwrite(strcat(file_dir,'\Unique-N-vdsl.txt'), frame, 'newline','pc','-append')
                        end
                    end
                end
            %end
        end

    elseif strcmp('Tp',var_name)    %unique in Tp values
        for Tp=1:64                 %possible Tp vals
            if Tp_uv(Tp)==0
               % for i=1:5           %pointer to Mp vals
               i=max(1,min(5,round((4*length(find(Tp_uv(1:64)))+63)/63)));     
               Mp=Mp_all(i);
                    perb=((Tp*n)/Mp)*floor((q1*Mp)/(Tp*n)); %bytes
                    per= (8*perb*1000)/(Lp*4000);
                    if per<=20 && per>=1 && mod(Tp,Mp)==0   %per lies between 1ms and 20ms
                        Gp=ceil(((msgp*Tp*s*1000)/(8*Mp*4000))+((6*n*Tp)/(perb*Mp)));
                        if Gp<=32 && Gp>=1
                            if (inp_co_uv(frame(1))<inp_co(frame(1)))
                                if Tp_uv(Tp)==0                 %unique occurance of Tp?
                                    if (Mp/s)<=64               %rule1
                                        de=(ceil(Mp/s));
                                        ce=mod(Gp,Tp);          %rule2
                                        if (floor(Gp/Tp)*de)+(floor(de/Tp)*(ce))+min(mod(de,Tp),ce)<=8 
                                            if mod((n-r),Mp)==0 %Bpn-Nfec integer relation
                                                chk=1;
                                                op(1:32)=0;
                                                for x=1:Tp
                                                    if x<=((Gp)-(Tp*floor(Gp/Tp)))
                                                        op(x)=ceil(Gp/Tp);
                                                    else op(x)=floor(Gp/Tp);
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end

                                if chk==1
                                    chk=0;
                                    Gp_uv(Gp)=1;    %flags
                                    Tp_uv(Tp)=1;
                                    D_uv(d)=1;
                                    N_uv(n)=1;
                                    Mp_uv(i)=1;
                                    inp_co_uv(frame(1))=inp_co_uv(frame(1))+1; %increase INP counts

                                    frame(8)=msgp;
                                    frame(9)=Tp;
                                    frame(10)=Mp;
                                    frame(11)=Gp;
                                    frame(12)=perb;
                                    frame(13)=per;
                                    frame(14)=(perb*Mp*Gp)/(n*Tp);  %SEQ
                                    frame(15:46)=op(1:32);

                                    dlmwrite(strcat(file_dir,'\Unique-Tp-vdsl.txt'), frame, 'newline','pc','-append')
                                end
                            else return     
                            end
                        end
                    end
               % end
            end
        end
    elseif strcmp('Gp',var_name)
        for Tp=1:64             %possible Tp vals
            for i=1:5           %pointer to Mp vals
                Mp=Mp_all(i);
                perb=((Tp*n)/Mp)*floor((q1*Mp)/(Tp*n)); %bytes
                per= (8*perb*1000)/(Lp*4000);
                if per<=20 && per>=1 && mod(Tp,Mp)==0   %per lies between 1ms and 20ms
                    Gp=ceil(((msgp*Tp*s*1000)/(8*Mp*4000))+((6*n*Tp)/(perb*Mp)));
                    if Gp<=32 && Gp>=1 && Gp_uv(Gp)==0
                        if (inp_co_uv(frame(1))<inp_co(frame(1)))
                            if  Gp_uv(Gp)==0                %unique Gp value?
                                if (Mp/s)<=64               %rule1
                                    de=(ceil(Mp/s));
                                    ce=mod(Gp,Tp);          %rule2
                                    if (floor(Gp/Tp)*de)+(floor(de/Tp)*(ce))+min(mod(de,Tp),ce)<=8 
                                        if mod((n-r),Mp)==0 %Bpn-Nfec integer relation
                                            chk=1;op(1:32)=0;
                                            for x=1:Tp
                                                if x<=((Gp)-(Tp*floor(Gp/Tp)))
                                                    op(x)=ceil(Gp/Tp);
                                                else op(x)=floor(Gp/Tp);
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                            if chk==1
                                chk=0;
                                Gp_uv(Gp)=1;    %flags
                                Tp_uv(Tp)=1;
                                D_uv(d)=1;
                                N_uv(n)=1;
                                Mp_uv(i)=1;
                                inp_co_uv(frame(1))=inp_co_uv(frame(1))+1;  %increment count of INP occurance

                                frame(8)=msgp;
                                frame(9)=Tp;
                                frame(10)=Mp;
                                frame(11)=Gp;
                                frame(12)=perb;
                                frame(13)=per;
                                frame(14)=(perb*Mp*Gp)/(n*Tp);%SEQ
                                frame(15:46)=op(1:32);

                                dlmwrite(strcat(file_dir,'\Unique-Gp-vdsl.txt'), frame, 'newline','pc','-append')
                            end
                            
                        else return
                        end
                    end
                end
            end
        end
    elseif strcmp('Mp',var_name)
        for Tp=1:64             %possible Tp vals
            for i=1:5           %pointer to Mp vals
                Mp=Mp_all(i);
                if Mp_uv(i)==0
                    perb=((Tp*n)/Mp)*floor((q1*Mp)/(Tp*n)); %bytes
                    per= (8*perb*1000)/(Lp*4000);
                    if per<=20 && per>=1 && mod(Tp,Mp)==0   %per lies between 1ms and 20ms
                        Gp=ceil(((msgp*Tp*s*1000)/(8*Mp*4000))+((6*n*Tp)/(perb*Mp)));
                        if  Mp_uv(i)==0                     %unique occurance of Mp?
                            if Gp<=32 && Gp>=1
                                if (inp_co_uv(frame(1))<inp_co(frame(1)))
                                    if (Mp/s)<=64               %rule1
                                        de=(ceil(Mp/s));
                                        ce=mod(Gp,Tp);          %rule2
                                        if (floor(Gp/Tp)*de)+(floor(de/Tp)*(ce))+min(mod(de,Tp),ce)<=8 
                                            if mod((n-r),Mp)==0 %bpn-Nfec integer condition
                                                chk=1;op(1:32)=0;
                                                for x=1:Tp
                                                    if x<=((Gp)-(Tp*floor(Gp/Tp)))
                                                        op(x)=ceil(Gp/Tp);
                                                    else op(x)=floor(Gp/Tp);
                                                    end
                                                end
                                            end
                                        end
                                    end
                                    if chk==1
                                        chk=0;
                                        Gp_uv(Gp)=1;
                                        Tp_uv(Tp)=1;
                                        D_uv(d)=1;
                                        N_uv(n)=1;
                                        Mp_uv(i)=1;
                                        inp_co_uv(frame(1))=inp_co_uv(frame(1))+1;

                                        frame(8)=msgp;
                                        frame(9)=Tp;
                                        frame(10)=Mp;
                                        frame(11)=Gp;
                                        frame(12)=perb;
                                        frame(13)=per;
                                        frame(14)=(perb*Mp*Gp)/(n*Tp);%SEQ
                                        frame(15:46)=op(1:32);

                                        dlmwrite(strcat(file_dir,'\Unique-Mp-vdsl.txt'), frame, 'newline','pc','-append')

                                    end
                                else return
                                end
                            end
                            
                        end
                    end
                end
            end
        end
    end
elseif strcmp('adsl',sl)
    if strcmp('D',var_name) %unique D values, adsl case
        for Tp=1:64         %Table 7-8/G.992.3
            for i=1:5
                Mp=Mp_all(i);
                if (s>=(Mp/16))&& (s<=(32*Mp))          %Sp-Mp conditon, Table 7-8/G.992.3
                    OR=(Mp*32*a)/(Tp);                  %kbit/s, Table 7-7/G.992.3
                    if OR>=0.8 && OR<=64 && OR >msgp    %Overhead Rate Constraints, Table 7-8/G.992.3
                        msgc=ceil((6*msgp)/((OR)-msgp));%msgp=MSGmin, using rel from Framing_VDSL_V05f(version1).xls adsl2
                        SEQ=msgc+6;                     %Table 7-7/G.992.3
                        per= (Tp*SEQ)/(4*a*Mp);         %Table 7-7/G.992.3
                        if per<=20 && per>=15           %Table 7-8/G.992.3,Overhead Channel Period
                            if D_uv(find(D_adsl==d))==0 %check if this value of D is unique?
                                if mod((n-r),Mp)==0 %ensuring Bpn-nfec interger condition, Figure 7-7/G.992.3
                                    Tp_uv(Tp)=1;                %flag to Tp
                                    D_uv(find(D_adsl==d))=1;    %flag to the corresponding D
                                    N_uv(n)=1;                  %flag to N
                                    Mp_uv(i)=1;                 %flag to Mp
                                    inp_co_uv(frame(1))=inp_co_uv(frame(1))+1;  %increase INP occurance count

                                    frame(8)=msgp;
                                    frame(9)=Tp;
                                    frame(10)=Mp;
                                    frame(11)=per;
                                    frame(12)=SEQ;
                                    frame(13)=OR;

                                    dlmwrite(strcat(file_dir,'\Unique-D-adsl.txt'), frame, 'newline','pc','-append')
                                end
                            else return
                            end
                        end
                    end
                end
            end
        end
    elseif strcmp('N',var_name) %unique N value,adsl
        for Tp=1:64
            %for i=1:5
                i=max(1,min(5,round((4*length(find(N_uv(32:255)))+223)/223)));
                Mp=Mp_all(i);
                if (s>=(Mp/16))&& (s<=(32*Mp))
                    OR=(Mp*32*a)/(Tp);
                    if OR>=0.8 && OR<=64 && OR >msgp
                        msgc=ceil((6*msgp)/((OR)-msgp));
                        SEQ=msgc+6;
                        per= (Tp*SEQ)/(4*a*Mp);
                        if per<=20 && per>=15
                            if N_uv(n)==0                       %Is this an unique occurance of N?
                                if mod((n-r),Mp)==0             %Bpn-Nfec integer condition 
                                    Tp_uv(Tp)=1;                %flags
                                    D_uv(find(D_adsl==d))=1;
                                    N_uv(n)=1;
                                    Mp_uv(i)=1;
                                    inp_co_uv(frame(1))=inp_co_uv(frame(1))+1;  %increase INP count

                                    frame(8)=msgp;
                                    frame(9)=Tp;
                                    frame(10)=Mp;
                                    frame(11)=per;
                                    frame(12)=SEQ;
                                    frame(13)=OR;

                                    dlmwrite(strcat(file_dir,'\Unique-N-adsl.txt'), frame, 'newline','pc','-append')
                                end
                            else return
                            end
                        end
                    end
                end
            %end
        end
    elseif strcmp('Tp',var_name)                                %unique Tp value,adsl
        for Tp=1:64
            if Tp_uv(Tp)==0                                     %is this unique occurance of Tp
                %for i=1:5
                    i=max(1,min(5,round((4*length(find(Tp_uv(1:64)))+63)/63))); 
                    Mp=Mp_all(i);
                    if (s>=(Mp/16))&& (s<=(32*Mp))
                        OR=(Mp*32*a)/(Tp);
                        if OR>=0.8 && OR<=64 && OR >msgp
                            msgc=ceil((6*msgp)/((OR)-msgp));
                            SEQ=msgc+6;
                            per= (Tp*SEQ)/(4*a*Mp);
                            if per<=20 && per>=15
                                if (inp_co_uv(frame(1))<inp_co(frame(1)))
                                    if mod((n-r),Mp)==0             %bpn-Nfec condition
                                        Tp_uv(Tp)=1;                %flags
                                        D_uv(find(D_adsl==d))=1;
                                        N_uv(n)=1;
                                        Mp_uv(i)=1;
                                        inp_co_uv(frame(1))=inp_co_uv(frame(1))+1;  %increase INP occurance counts

                                        frame(8)=msgp;
                                        frame(9)=Tp;
                                        frame(10)=Mp;
                                        frame(11)=per;
                                        frame(12)=SEQ;
                                        frame(13)=OR;

                                        dlmwrite(strcat(file_dir,'\Unique-Tp-adsl.txt'), frame, 'newline','pc','-append')
                                    end
                                else return
                                end
                            end
                        end
                    end
                %end
            end
        end
    elseif strcmp('Mp',var_name)                                    %unique Mp value,adsl
        for Tp=1:64
            for i=1:5
                Mp=Mp_all(i);
                if (s>=(Mp/16))&& (s<=(32*Mp))
                    if Mp_uv(i)==0                                  %is it an unique occurance of Mp?
                        OR=(Mp*32*a)/(Tp);
                        if OR>=0.8 && OR<=64 && OR >msgp
                            msgc=ceil((6*msgp)/((OR)-msgp));
                            SEQ=msgc+6;
                            per= (Tp*SEQ)/(4*a*Mp);
                            if per<=20 && per>=15
                                if (inp_co_uv(frame(1))<inp_co(frame(1)))
                                    if mod((n-r),Mp)==0             %bpn-Nfec condition
                                        Tp_uv(Tp)=1;                %flags
                                        D_uv(find(D_adsl==d))=1;
                                        N_uv(n)=1;
                                        Mp_uv(i)=1;
                                        inp_co_uv(frame(1))=inp_co_uv(frame(1))+1;  %increase INP occurance counts

                                        frame(8)=msgp;
                                        frame(9)=Tp;
                                        frame(10)=Mp;
                                        frame(11)=per;
                                        frame(12)=SEQ;
                                        frame(13)=OR;

                                        dlmwrite(strcat(file_dir,'\Unique-Mp-adsl.txt'), frame, 'newline','pc','-append')
                                    end
                                else return
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end