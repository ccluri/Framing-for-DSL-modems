function INPrate0_gui(file_name)
% This generates all the possible INPrate combinations for a given min
% INP rate and max delay over a range of L, for different freq plans of vdsl
% and adsl and adsl+ limits for the same are assumed to be available from
% the gui_inp.m file (the values entered in the gui)
%
% An addition possibility of evaluating INPrates when at fast path ie.,
% when r=0, d=1 case.
%
% the default user defined limits are part of the global variable named
% looplimits. the start,stop,step sizes of Min INP, Max delay, R and Q(in
% case of vdsl), max, min,step sizes of L, smax,smin, min_msg are part
% of this, in addition of some checks for fast path,framing, and choice of
% frequency plans

% Algorithm of the INPrate.
% the limits of nfec nmax and nmin of the particular L are set based on
% the limits of s, corresponding to the freq plan. (if they exceed the
% respective min-max limits of nfec, 32-255, are chosen)
% limiting case: nfec=s*L/8; nmax:smax, nmin:smin ; & nfec range= [32,255];
%
% All the possible cofigarations for R,q combinations are determied and
% saved in the variable iq_para. no.of possible combinations = num_inp_sets
% limiting case: floor(r/(2*q))=0; d=(inpmin*q)/(floor(r/2q)*8);
% divide by zero error.
%
% evaluate d=(inp_min*q)/(floor(r/2q)*8)
% if d=1 (force make d=2)
% limiting case: d>dmax of freq plan; mark as invalid configation
%
% evaluate I(=J)=floor(Fs*delmax*L/8); & nfec=q*J; (for valid cases)
% check for limiting cases (in this sequence)
% 1)if nfec>nmax, force J to limiting case of nmax, nfec=nmax
% 2)if (d-1)*(J-1)>max_interlvr, force J to maxinterlever case, eval nfec
% 3)if nfec<nmin; mark as invalid configaration
%
% finding the biggest interlever depth possible (for valid cases)
% while d<dmax
%   checking for coprime condition of (d,n),update if through.
%   increment d and check for limiting cases of
%   1)evaluate J, nfec_new 2)if nfec_new>nmax, force J by nfec_new=nmax
%   3)if (J-1)(d-1)>max interlvr, force J to maxinterlvr case,eval nfec_new
%   4)if nfec_new <nmin, mark as invalid configaration
%   5)if nfec_new <nfec, update nfec=nfec_new,and continue(try for nxt d)
%
% output generation
% evaluate: actual (s, inp, delay) tdr and adr.
% save all valid cases iq_max into the variable named output.
% if needed to be framed, pass the output to frame_gui funtion. otherwise
% print into the file_name.txt opened (passed from gui_inp)


global looplimits %input values for the limits

fid = fopen(file_name, 'wt'); %open the file to put the resulting output

N_R_Q_PARA      = 36;   % number of possible R, q combinations
VALID           = 0;	% result in inp_set struct is valid
INVALID         = -1;	% result in inp_set struct is not valid

% LIMIT Codes iq_para(12)
LIMIT_NONE      = 0;
D_MAX           = 1;
D_MIN           = 1;
S_MAX           = 16;
S_MIN           = 32;
FEC_MAX         = 256;
FEC_MIN         = 512;
INTL_MEM_MAX    = 1024;
AGG_MAX         = 4096;
NO_COPRIME      = 32768;

% Reed-Solomon check bytes
R_MIN           =	2;	% smallest number of check bytes - must be an even value
R_MAX           =	16;	% highest number of check bytes  - must be an even value

% GCI (generalized convolutional interleaver parameter)
Q_MIN           =	1;	% minimum value for q parameter (Nfec = q * I; max. range of q: 1 <=q <= 8)
Q_MAX           =	8;	% maximum value for q parameter (Nfec = q * I; max. range of q: 1 <=q <= 8)

% Code word size range
NFEC_MIN        = 32;
NFEC_MAX        = 255;

%% Set and Check the config_data

INP_start   = looplimits(1);        %user defined
INP_end     = looplimits(2);
INP_step    = looplimits(3);


DEL_start   = looplimits(4)/1000;   %user defined
DEL_end     = looplimits(5)/1000;   %converting into seconds
DEL_step    = looplimits(6)/1000;

L_min       = looplimits(7);        %user defined
L_max       = looplimits(8);
L_step      = looplimits(9);

R_min       = looplimits(10);       %user defined
R_max       = looplimits(11);

Q_min       = looplimits(12);       %user defined
Q_max       = looplimits(13);

D_max       = looplimits(14);       %profile defined


inv_S_min   = looplimits(15);       %profile defined
S_min       = 1.0/inv_S_min;
S_max       = looplimits(16);

IntlMem_max = looplimits(17);       %profile defined

Fs          = looplimits(18);       %profile defined
msg_min     = looplimits(19);       %user defined


if looplimits(20)==1        %adsl/vdsl
    sl='vdsl';
elseif looplimits(20)==2
    sl='adsl';
end

if looplimits(21)==1        %fast path case
    intlv='f';
elseif looplimits(21)==2    %full path case
    intlv='o';
end

% if looplimits(22)==1  %stream
%     stream='ds';
% elseif looplimits(22)==2
%     stream='us';
% end

if looplimits(23)==1
    frame=1;            %frame
else frame=0;           %do not frame
end

%% Configure profile limits
limits = struct('Smax',S_max, 'Smin',S_min, 'Dmax',D_max, 'IntlMem',IntlMem_max, 'NFECmax',NFEC_MAX, 'NFECmin',NFEC_MIN);
limits.IntlMem;
iq_max = zeros(13,1);
iq_para = zeros (13,N_R_Q_PARA);

% the members are as follows iq_para
%1.  R          long    number of checkbytes
%2.  q          long    q parameter of GCI
%3.  R_2q       long    floor(R/2q)
%4.  L          long    bit per symbol
%5.  D          long    resulting interleaver depth of GCI
%6.  Nfec       long    resulting code word size
%7.  Nfec_min   long    min code word size (limits on Smin, L)
%8.  Nfec_max   long    max code word size (limits on Smax, L)
%9.  inp        double  actual inp
%10. del        double  actual delay
%11. adr        double  actual aggregate data rate
%12. limit      long    optimized parameters constrained by profile limits
%13. stat       long    parameter set status: 0: valid, <> 0: not valid

if intlv~='f'   %full path

    for INPmin = INP_start:INP_step:INP_end

        for DELmax = DEL_start:DEL_step:DEL_end

            iq_max(11) = 0.0;
            iq_max(13) = INVALID;

            for L = L_min:L_step:L_max

                Nmax = floor(limits.Smax * L/8.0); %set limits to nfec
                if Nmax > limits.NFECmax
                    Nmax = limits.NFECmax;
                end

                Nmin = ceil(limits.Smin * L/ 8.0); %set limits to nfec
                if Nmin < limits.NFECmin
                    Nmin = limits.NFECmin;
                end

                num_inp_sets = 0;

                for r = R_max:-2:R_min
                    for q = Q_min:1:min(r/2,Q_max)                      %ensuring floor(r/2q)!=0

                        iq_para(num_inp_sets*13+1)  = r;
                        iq_para(num_inp_sets*13+2)  = q;
                        iq_para(num_inp_sets*13+3)  = floor(r/(2*q));
                        iq_para(num_inp_sets*13+4)  = L;
                        iq_para(num_inp_sets*13+5)  = 0;    %d
                        iq_para(num_inp_sets*13+6)  = 0;    %nfec
                        iq_para(num_inp_sets*13+7)  = Nmin;
                        iq_para(num_inp_sets*13+8)  = Nmax;
                        iq_para(num_inp_sets*13+9)  = 0.0;  % actual INP
                        iq_para(num_inp_sets*13+10) = 0.0;  % actual delay
                        iq_para(num_inp_sets*13+11) = 0.0;  % actual ADR
                        iq_para(num_inp_sets*13+12) = LIMIT_NONE;
                        iq_para(num_inp_sets*13+13) = VALID;
                        num_inp_sets = num_inp_sets + 1;

                    end
                end


                for i = 0:1:num_inp_sets-1

                    d = ceil(INPmin * iq_para(i*13+4) / iq_para(i*13+3) / 8.0); %ceil(INPmin*L/(8*flr(r/2q)))

                    if sl=='adsl'
                        d=getD_adsl(d,looplimits(22));
                    end

                    if d <= limits.Dmax

                        if d == 1
                            d = 2;	% force use of interleaver
                            if mod(floor(iq_para(i*13+12)/D_MIN), 2) == 0
                                iq_para(i*13+12) = iq_para(i*13+12) + D_MIN;
                            end
                        end

                        iq_para(i*13+5) = d;


                    else            % required interleaver depth is out of range
                        iq_para(i*13+13) = INVALID;
                        if mod(floor(iq_para(i*13+12)/D_MAX), 2) == 0
                            iq_para(i*13+12) = iq_para(i*13+12) + D_MAX;
                        end
                    end

                end


                for i = 0:1:num_inp_sets-1

                    if iq_para(i*13+13) == INVALID
                        continue;
                    end

                    J = floor (Fs/8.0 * DELmax * iq_para(i*13+4) / (iq_para(i*13+5)-1)  + 1.0); %I=floor(8Fs*del*L/(D-1)+1)
                    nfec = J * iq_para(i*13+2);

                    if nfec > iq_para(i*13+8)       % check whether code word size exceeds upper limit	(255, S related limit)
                        J = floor(iq_para(i*13+8) / iq_para(i*13+2));       % new (smaller) J which forces nfec into upper limits
                        nfec = J * iq_para(i*13+2); % nfec =nmax
                        if mod(floor(iq_para(i*13+12)/FEC_MAX), 2) == 0
                            iq_para(i*13+12) = iq_para(i*13+12) + FEC_MAX;	% report limit condition
                        end
                    end

                    if limits.IntlMem < ((iq_para(i*13+5) - 1)*(J - 1))             % check for overall interleaver memory limit
                        J = floor (limits.IntlMem / (iq_para(i*13+5) - 1) + 1.0);   % Select J which does not exceed interleaver mem
                        nfec = J * iq_para(i*13+2);                                 % N for Jmax
                        if mod(floor(iq_para(i*13+12)/INTL_MEM_MAX), 2) == 0
                            iq_para(i*13+12) = iq_para(i*13+12) + INTL_MEM_MAX;     % report limit condition
                        end
                    end


                    if (nfec < iq_para(i*13+7))     % if code word size nfec is smaller than Nfec_min --> R,q option not feasible.
                        iq_para(i*13+13) = INVALID; % mark result as invalid, still take result as Nfec
                        if mod(floor(iq_para(i*13+12)/FEC_MIN), 2) == 0
                            iq_para(i*13+12) = iq_para(i*13+12) + FEC_MIN;	% report limit condition
                        end
                    end

                    iq_para(i*13+6) = nfec;         % save code word size (even if it is out of range)

                end

                for i = 0:1:num_inp_sets-1

                    if iq_para(i*13+13) == INVALID
                        continue;
                    end

                    nfec = iq_para(i*13+6);         %I=nfec/q{=iq_para(i*13+2)};
                    d = iq_para(i*13+5);


                    while d <= limits.Dmax
                        % CheckCoPrime of d and I, found valid result
                        s_temp1 = d;
                        s_temp2 = nfec/iq_para(i*13+2); %=J
                        while s_temp2 > 0   % unless remainder = 0
                            s_temp3 = floor(s_temp1 / s_temp2);
                            s_temp4 = s_temp1 - s_temp2 * s_temp3;
                            s_temp1 = s_temp2;
                            s_temp2 = s_temp4;
                        end

                        if s_temp1 == 1     % co-prime
                            iq_para(i*13+6) = nfec;
                            iq_para(i*13+5) = d;
                            iq_para(i*13+9) =  (8.0 * d * iq_para(i*13+3)) / iq_para(i*13+4);                           %actual INP
                            iq_para(i*13+10) = 8.0/Fs * ((d - 1) / iq_para(i*13+4)) * (nfec / iq_para(i*13+2) - 1.0);   %actual delay
                            iq_para(i*13+11) = Fs * iq_para(i*13+4) * (1.0 - iq_para(i*13+1) / nfec);                   %actual adr
                            break;
                        end

                        d = d + 1; % try next bigger interleaver depth

                        if sl=='adsl'
                            d=getD_adsl(d,looplimits(22));
                        end

                        % compute code word size for new d. Options: a) still the same nfec b)smaller nfec
                        J = floor (( Fs/8.0 * DELmax * iq_para(i*13+4)) / (d - 1)  + 1.0);
                        nfec_new = J * iq_para(i*13+2); % initial code word size (may exceed the range for valid nfec at this point)

                        if nfec_new > iq_para(i*13+8)
                            J = floor(iq_para(i*13+8) / iq_para(i*13+2));       % new (smaller) J which forces nfec into upper limits with new d
                            nfec_new = J * iq_para(i*13+2);
                            if mod(floor(iq_para(i*13+12)/FEC_MAX), 2) == 0
                                iq_para(i*13+12) = iq_para(i*13+12) + FEC_MAX;	% report limit condition
                            end
                        end

                        if limits.IntlMem < ((d - 1)*(J - 1))           % check for overall interleaver memory limit
                            J = floor (limits.IntlMem / (d - 1) + 1.0); % Select J which does not exceed interleaver mem
                            nfec_new = J * iq_para(i*13+2);
                            if mod(floor(iq_para(i*13+12)/INTL_MEM_MAX), 2) == 0
                                iq_para(i*13+12) = iq_para(i*13+12) + INTL_MEM_MAX;	% report limit condition
                            end
                        end


                        if nfec_new < iq_para(i*13+7)   % if new code word size is below the limit of nfec --> R,q option not feasible.
                            nfec = nfec_new;
                            iq_para(i*13+13) = INVALID;	% mark result as invalid, still take result as Nfec
                            if mod(floor(iq_para(i*13+12)/FEC_MIN), 2) == 0
                                iq_para(i*13+12) = iq_para(i*13+12) + FEC_MIN;	% report limit condition
                            end
                            break;
                        end

                        if nfec_new < nfec	% if new cord word size is less than the previous one start again with initial d (tab[i].D)
                            nfec = nfec_new;
                            d = iq_para(i*13+5);
                        end

                    end

                    if d > limits.Dmax

                        iq_para(i*13+13) = INVALID;                             % mark result as invalid, still take result as Nfec
                        if mod(floor(iq_para(i*13+12)/NO_COPRIME), 2) == 0
                            iq_para(i*13+12) = iq_para(i*13+12) + NO_COPRIME;	% report limit condition
                        end
                    end
                end

                %% check if there is a successful configuration and if yes then save the optimal configuration
                % void evaluate_result (int num, struct inp_set tab[], struct inp_set opt_set[])
                % evaluate_result (num_inp_sets, iq_para, iq_max);

                i = 0; %pointer to iq_para

                while i < num_inp_sets

                    if iq_para(i*13+13) == VALID

                        iq_max = iq_para(:,i+1);

                        S= 8.0 * iq_max(6) / iq_max(4);	% compute S
                        tdr = iq_max(4) * Fs;           % compute Total Data Rate

                        output(1)=INPmin;
                        output(2)=DELmax;
                        output(3)=iq_max(4);            %Lp
                        output(4)=iq_max(1);            %r
                        output(5)=iq_max(2);            %q
                        output(6)=iq_max(5);            %d
                        output(7)=iq_max(6);            %nfec
                        output(8)=1.0/S;                %1/s
                        output(9)=iq_max(9);            %inp
                        output(10)=iq_max(10)*1000.0;   %del*1000
                        output(11)=tdr*1e-6;            %tdr
                        output(12)=iq_max(11)*1e-6;     %adr

                        if (frame)
                            frame_gui(Fs,sl,msg_min,output,fid); %Lp,1/s,Fs,r,nfec,adsl or vdsl,msgp_min
                        else fprintf(fid,'%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n',output(1:12));
                        end

                    end
                    i=i+1;
                end
            end
        end
    end
end


if intlv=='f' %fast path case

    %     r=0;
    %     q=1;
    %     d=1;
    nfec=255;

    for INPmin = INP_start:INP_step:INP_end
        for DELmax = DEL_start:DEL_step:DEL_end
            for L = L_min:L_step:L_max
                % for r=0:2:R_MAX
                for q=1:1:8

                    if r==0
                        d=1;
                        inp =  0;   % actual INP
                    else
                        d = ceil(INPmin * L/ floor(r/(2*q)) / 8.0);
                        inp =  (8.0 * d * floor(r/(2*q))) /L;   % actual INP
                    end

                    if d==1
                        S= 8.0 * nfec / L;	% compute S
                        tdr = L * Fs;       % compute Total Data Rate
                        adr = Fs * L * (1.0 / nfec);
                        del = 0;
                        output(1)=INPmin;
                        output(2)=DELmax;
                        output(3)=L;
                        output(4)=r;
                        output(5)=q;
                        output(6)=d;
                        output(7)=nfec;
                        output(8)=1.0/S;
                        output(9)=inp;
                        output(10)=del*1000.0;
                        output(11)=tdr*1e-6;
                        output(12)=adr*1e-6;

                        if (frame)
                            frame_gui(Fs,sl,msg_min,output,fid); %Lp,1/s,Fs,r,nfec,adsl or vdsl,msgp_min
                        else fprintf(fid,'%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n',output(1:12));
                        end

                    end
                end
                %end
            end
        end
    end

end

fclose(fid);


clear looplimits

function [d]=getD_adsl(d,stream)    %adsl D selection, Table 7-8/G.992.5
if stream==1 %downstream
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





