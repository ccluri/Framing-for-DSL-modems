function Redundancy(stream)
%stream input is string input from of the generated cases, vdsl or adsl
%cases. This function is called after 'All' Unique cases from GUI are
%selected, The function removes redundacy in the down-up path.

if stream==1
    sl='-vdsl.txt';
else sl='-adsl.txt';
end

if exist (strcat('Unique-Mp',sl),'file');   %check if file exists
    m=dlmread(strcat('Unique-Mp',sl));      %read contents
    mf=1;                                   %check its existence
else mf=0; m=[];
end
if exist (strcat('Unique-Gp',sl),'file');
    g=dlmread(strcat('Unique-Gp',sl));
    gf=1;
else gf=0; g=[];
end
if exist (strcat('Unique-Tp',sl),'file');
    t=dlmread(strcat('Unique-Tp',sl));
    tf=1;
else tf=0; t=[];
end
if exist (strcat('Unique-N',sl),'file');
    n=dlmread(strcat('Unique-N',sl));
    nf=1;
else nf=0; n=[];
end
if exist (strcat('Unique-D',sl),'file');
    d=dlmread(strcat('Unique-D',sl));
    df=1;
else df=0; d=[];
end

r_all=[m;g;t;n;d];  %all possible unique values before removing redundancy
all=[];

if (mf)
    for i=1:size(m,1)
        if(gf)
            [x,y]=find(g(:,11)==m(i,11));   %check for redundancy in Gp output (already covered in Mp)
            if (~isempty(x))                %if redundancy exists
                g(x,:)=[];                  %remove redundancy
                r_all=[];
                r_all=[m;g;t;n;d];          %update unique cases
            end
        end
        if (tf)
            [x,y]=find(t(:,9)==m(i,9));
            if (~isempty(x))                    %check for redundancy in Tp output (already covered in Mp)
                if length(find(r_all(:,11)==t(x,11)))>1 %ensure occurance of this Tp has not a unique Gp
                    t(x,:)=[];                  %remove redundancy
                    r_all=[];
                    r_all=[m;g;t;n;d];          %update unique cases
                end
            end
        end
        if (nf)
            [x,y]=find(n(:,3)==m(i,3));
            if (~isempty(x))                    %check for redundancy in N output (already covered in Mp)
                if length(find(r_all(:,11)==n(x,11)))>1 && length(find(r_all(:,9)==n(x,9)))>1   %ensure occurance of this N is not a unique Tp or Gp
                    n(x,:)=[];                  %remove redundancy
                    r_all=[];
                    r_all=[m;g;t;n;d];          %update unique cases
                end
            end
        end
        if (df)
            [x,y]=find(d(:,5)==m(i,5));
            if (~isempty(x))                    %check for redundancy in D output (already covered in Mp)
                if length(find(r_all(:,11)==d(x,11)))>1 && length(find(r_all(:,9)==d(x,9)))>1 && length(find(r_all(:,3)==d(x,3)))>1
                    %ensure occurance of this D is not a unique Tp or Gp or N
                    d(x,:)=[];                  %remove redundancy
                    r_all=[];
                    r_all=[m;g;t;n;d];          %update unique cases
                end
            end
        end
    end

end

if (gf)
    for i=1:size(g,1)
        if (tf)                                 %check for redundancy in Tp output (already covered in Gp)
            [x,y]=find(t(:,9)==g(i,9));
            if (~isempty(x))                    %ensure occurance of this Tp has not a unique Gp
                if length(find(r_all(:,11)==t(x,11)))>1
                    t(x,:)=[];                  %remove redundancy
                    r_all=[];
                    r_all=[m;g;t;n;d];          %update unique cases
                end
            end
        end
        if (nf)
            [x,y]=find(n(:,3)==g(i,3),1,'first');
            if (~isempty(x))                    %check for redundancy in N output (already covered in Gp)
                if length(find(r_all(:,11)==n(x,11)))>1 && length(find(r_all(:,9)==n(x,9)))>1 %ensure occurance of this Tp has not a unique Gp
                    n(x,:)=[];                   %remove redundancy
                    r_all=[];
                    r_all=[m;g;t;n;d];           %update unique cases
                end
            end
        end
        if (df)
            [x,y]=find(d(:,5)==g(i,5));
            if (~isempty(x))                    %check for redundancy in D output (already covered in Gp)
                if length(find(r_all(:,11)==d(x,11)))>1 && length(find(r_all(:,9)==d(x,9)))>1 && length(find(r_all(:,3)==d(x,3)))>1
                    %ensure occurance of this D is not a unique Tp or Gp or N
                    d(x,:)=[];                  %remove redundancy
                    r_all=[];
                    r_all=[m;g;t;n;d];          %update unique cases
                end
            end
        end
    end

end

if (tf)
    for i=1:size(t,1)
        if (nf)                                 %check for redundancy in N output (already covered in Tp)
            [x,y]=find(n(:,3)==t(i,3),1,'first');
            if (~isempty(x))
                if length(find(r_all(:,11)==n(x,11)))>1 && length(find(r_all(:,9)==n(x,9)))>1 %ensure occurance of this Tp has not a unique Gp
                    n(x,:)=[];                  %remove redundancy
                    r_all=[];
                    r_all=[m;g;t;n;d];          %update unique cases
                end
            end
        end
        if (df)                                 %check for redundancy in D output (already covered in Tp)
            [x,y]=find(d(:,5)==t(i,5));
            if (~isempty(x))
                if length(find(r_all(:,11)==d(x,11)))>1 && length(find(r_all(:,9)==d(x,9)))>1 && length(find(r_all(:,3)==d(x,3)))>1
                    %ensure occurance of this D is not a unique Tp or Gp or N
                    d(x,:)=[];                  %remove redundancy
                    r_all=[];
                    r_all=[m;g;t;n;d];          %update unique cases
                end
            end
        end
    end

end

if (nf)                                         %check for redundancy in D output (already covered in N)
    for i=1:size(n,1)
        if (df)
            [x,y]=find(d(:,5)==n(i,5),1,'first');
            if (~isempty(x))
                if length(find(r_all(:,11)==d(x,11)))>1 && length(find(r_all(:,9)==d(x,9)))>1 && length(find(r_all(:,3)==d(x,3)))>1
                    %ensure occurance of this D is not a unique Tp or Gp or N
                    d(x,:)=[];                  %remove redundancy
                    r_all=[];
                    r_all=[m;g;t;n;d];          %update unique cases
                end
            end
        end
    end

end

all=r_all;

if stream==1                    %formatting and printing into excel sheet.
    n_all(:,1:2)=all(:,1:2);    %inp,Lp
    n_all(:,3)=all(:,10);       %mp
    n_all(:,4)=all(:,9);        %Tp
    n_all(:,5)=all(:,8);        %min_msg
    n_all(:,6)=floor((all(:,3)-all(:,7))./all(:,10))-ceil(all(:,11)./all(:,9));%Bp=floor((n-r)/m)-ceil(gp/tp)
    n_all(:,7)=all(:,13);       %per
    n_all(:,8)=all(:,14);       %seq
    n_all(:,9)=all(:,5);        %d
    n_all(:,10)=all(:,6);       %q
    n_all(:,11)=all(:,3)./all(:,6);%I=n/q
    n_all(:,12)=all(:,7);       %R
    n_all(:,13)=all(:,3);       %N
    n_all(:,14)=all(:,4);       %s
    n_all(:,15)=all(:,11);      %Gp
    
    n_all(:,16)=4*n_all(:,2);   %TDR =Fs*Lp
    n_all(:,17)=4*8*n_all(:,15).*n_all(:,3)./(n_all(:,14).*n_all(:,4));     %OR = 8*Fs*Gp*Mp/(s*Tp)
    n_all(:,18)=n_all(:,16).*(1-(n_all(:,12)./n_all(:,13)))-n_all(:,17);     %NDR = TDR*(1-(R/Nfec))-OR
    n_all(:,19)=(n_all(:,9)-1).*(n_all(:,11)-1);    %Actual_Intlv
   
    n_all(:,20:51)=all(:,15:46);%opi
    
    names={'INP','Lp','Mp','Tp','msg','Bp','per','seq','D','q','I','R','N','s','Gp','TDR','OR','NDR','Act_Intlv','Opi'};
    xlswrite('vdsl_uniquedata', names,'Sheet1','A1');
    xlswrite('vdsl_uniquedata', n_all,'Sheet1','A2');
else                            %adsl case
    n_all(:,1:2)=all(:,1:2);    %inp,Lp
    n_all(:,3)=all(:,10);       %mp
    n_all(:,4)=all(:,9);        %Tp
    n_all(:,5)=all(:,8);        %min_msg
    n_all(:,6)=floor((all(:,3)-all(:,7))./all(:,10))-1; %bp=floor((n-r)/mp)-1
    n_all(:,7)=all(:,11);       %per
    n_all(:,8)=all(:,12);       %seq
    n_all(:,9)=all(:,5);        %d
    n_all(:,10)=all(:,7);       %r
    n_all(:,11)=all(:,3);       %n
    n_all(:,12)=all(:,4);       %s
        
    n_all(:,13)=4*n_all(:,2);   %TDR 
    n_all(:,14)=8*4*n_all(:,3)./(n_all(:,4).*n_all(:,12)); %OR
    n_all(:,15)=((n_all(:,4).*(n_all(:,6)+1)-1).*n_all(:,3).*n_all(:,2)*4)./(n_all(:,4).*((n_all(:,3).*(n_all(:,6)+1))+n_all(:,10))); %NDR
    n_all(:,16)= (n_all(:,11)-1).*(n_all(:,9)-1); %Act_Intlv
    
    names={'INP','Lp','Mp','Tp','msg','Bp','per','seq','D','R','N','s','TDR','OR','NDR','Act_Intlv'};
    xlswrite('adsl_uniquedata', names,'Sheet1','A1');
    xlswrite('adsl_uniquedata', n_all,'Sheet1','A2');
end
