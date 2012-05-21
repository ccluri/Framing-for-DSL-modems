%PAIRING of DS and US unique cases of VDSL 

%input files to this script are Downstream and Upstream unique case excel
%sheets, The first prompt demands the downstream file followed by upstream
%the format of these files is (only for VDSL)
%***INP-Lp-Mp-Tp-msg-Bp-per-seq-D-q-I-R-N-s-Gp-TDR-OR-NDR-Act_Intlv-Opi*** 

%output is an excel sheet Master.xls, same format as Input files.(paired) 

%script gives priority to Upstream case. (all cases of upstream are met,
%with if necessary repetition of downstream cases) to change this 
%uncomment line 56 to 111 and comment 112-163, (downstream priority)

%plan specifications
Max_Intlv=65536  

%other specification of NDR
% max DS data rate 125
% max US data rate  75
% max DS+US data rate 200

disp('Select consolidated Downstream Unique case file')
[DS_filename, DS_pathname] = uigetfile({'*.xls'},'File Selector');
disp('Select consolidated Upstream Unique case file')
[US_filename, US_pathname] = uigetfile({'*.xls'},'File Selector');

DS=xlsread(strcat(DS_pathname,DS_filename));  
US=xlsread(strcat(US_pathname,US_filename));   

%colums in DS and US
%1:20 relevant fields. Remove 21:end Opi(i=2:32)bits
%change this to 2:21; 1 for profile ID
DS(:,21:end)=[];DS(:,2:21)=DS(:,1:20);DS(:,1)=0;
US(:,21:end)=[];US(:,2:21)=US(:,1:20);US(:,1)=0;

%save Actual Interleaver used (21st colum) in ascending order
DS=sortrows(DS,20);
US=sortrows(US,20);

n_ds=length(DS);
n_us=length(US);

%profile ID's after sorting
for i=1:n_ds
    DS(i,1)=i;
end
for i=1:n_us
    US(i,1)=i;
end

%defaulting zero value for indexes of flags
US(:,22:26)=0;
DS(:,22:26)=0;

DS(:,17:19)=DS(:,17:19)./1000;%Mbits/sec
US(:,17:19)=US(:,17:19)./1000;

names={'DS_ProfileID','INP','Lp','Mp','Tp','msg','Bp','per','seq','D','q','I','R','N',...
    's','Gp','TDR','OR','NDR','Act_Intlv','Opi','US_ProfileID','INP','Lp','Mp',...
    'Tp','msg','Bp','per','seq','D','q','I','R','N','s','Gp','TDR','OR',...
    'NDR','Act_Intlv','Opi'};

%starting values
Intlv=0;q=0;
NDR=0;
Intlv1=0;
NDR1=0;AM=0;AM1=0;
Max_flag_cnt=1;

i=1;
for j=n_us:-1:1
    if i==(n_ds+1)
        break;
    end
    if DS(i,20)+US(j,20)<=Max_Intlv
        if DS(i,19)<=125 && US(j,19)<=75 && (DS(i,19)+US(j,19))<=200
            if US(j,22)< Max_flag_cnt && DS(i,22)<Max_flag_cnt
                US(j,22)=US(j,22)+1;
                DS(i,22)=DS(i,22)+1;

                US(j,23)=i; %US number
                DS(i,24)=j; %DS number

                DS_US(i,1:26)=DS(i,:);
                DS_US(i,27:52)=US(j,:);

                i=i+1;
                if i==(n_ds+1)
                    break;
                end
            else AM=AM+1;
            end
        else  NDR=NDR+1;
        end
    else Intlv=Intlv+1;
    end
end

for Max_flag_cnt=2:100
    i=length(DS_US)+1;
    if i==(n_ds+1)
        break;
    end
        for j=n_us:-1:1
            if i<=length(DS)
            if DS(i,20)+US(j,20)<=Max_Intlv
                if DS(i,19)<=125 && US(j,19)<=75 && (DS(i,19)+US(j,19))<=200
                    if US(j,22)< Max_flag_cnt && DS(i,22)<Max_flag_cnt
                        US(j,22)=US(j,22)+1;
                        DS(i,22)=DS(i,22)+1;

                        US(j,23)=i; %US number
                        DS(i,24)=j; %DS number

                        DS_US(i,1:26)=DS(i,:);
                        DS_US(i,27:52)=US(j,:);

                        i=i+1;
                        if i==(n_ds+1)
                            break;
                        end
                    else AM=AM+1;
                    end
                else  NDR=NDR+1;
                end
            else Intlv=Intlv+1;
            end
        end
    end
end
output_final=horzcat(DS_US(:,1:21),DS_US(:,27:47)); 

% plot(DS_US(:,20),'.','color','r')
% hold on
% plot(DS_US(:,46),'.','color','b')
% plot(DS_US(:,46)+DS_US(:,20),'.','color','g')
% plot(1:2500,98304)
% hold off

% j=1;
% for i=n_ds:-1:1
%     if DS(i,19)<=125 && US(j,19)<=75 && (DS(i,19)+US(j,19))<=200
%         if DS(i,20)+US(j,20)<=Max_Intlv
%             if DS(i,22)< Max_flag_cnt && US(j,22)< Max_flag_cnt
%                 DS(i,22)=DS(i,22)+1;
%                 US(j,22)=US(j,22)+1;
% 
%                 US(j,25)=i; %US number
%                 DS(i,26)=j; %DS number
% 
%                 US_DS(j,1:26)=DS(i,:);
%                 US_DS(j,27:52)=US(j,:);
% 
%                 j=j+1;
%             end
%         else Intlv=Intlv+1;
%         end
%     else NDR=NDR+1;
%     end
% end
% 
% for Max_flag_cnt=2:100
%     j=length(US_DS)+1;
%     if j==(n_us+1)
%         break;
%     end
%     for i=n_ds:-1:1
%         if DS(i,19)<=125 && US(j,19)<=75 && (DS(i,19)+US(j,19))<=200
%             if DS(i,20)+US(j,20)<=Max_Intlv
%                 if DS(i,22)< Max_flag_cnt && US(j,22)< Max_flag_cnt
%                     DS(i,22)=DS(i,22)+1;
%                     US(j,22)=US(j,22)+1;
% 
%                     US(j,25)=i; %US number
%                     DS(i,26)=j; %DS number
% 
%                     US_DS(j,1:26)=DS(i,:);
%                     US_DS(j,27:52)=US(j,:);
% 
%                     j=j+1;
%                     if j==(n_us+1)
%                         break;
%                     end
%                 end
%             else Intlv=Intlv+1;
%             end
%         else NDR=NDR+1;
%         end
%     end
% end
% output_final=horzcat(US_DS(:,1:21),US_DS(:,27:47)); 
% 
% plot(US_DS(:,20),'.','color','r')
% hold on
% plot(US_DS(:,46),'.','color','b')
% plot(US_DS(:,46)+US_DS(:,20),'.','color','g')
% hold off

xlswrite('Master', names,'Sheet1','A1');
xlswrite('Master', output_final,'Sheet1','A2');
