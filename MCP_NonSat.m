
% How to run
% Example for N=10 and L=8:
% [THUAP, THM, THUS, PUAPloss, PMloss, PUSloss, DUAP, DM, DUS, JtM] = MCP_NonSat(10,8)
% all the DCF and MCP parameters are in this function: all_parameters.m


function [THUAP, THM, THUS, PUAPloss, PMloss, PUSloss, DUAP, DM, DUS, JtM] = MCP_NonSat(N,L)


    % US: Unicast senders
    % UAP: unicast of the AP
    % MAP: multicast of the AP

%     clc;
%     clear;
    %close all;
    %N=10;
    %L=7;


    %% Load all parameters
    [~, CW_min, sigma, DATA, ~, DIFS, ~, ~, delta, TUX, TUM, TCX, TCM, TXM, ~, ~, ~, ~, CW_max]=all_parameters();

    %% Simulation duration
    T = 1e9; % 1e11  default   T= 1e9 microsecond

    % arrival rate (packet/generic slot)
    lambdaUS = 10;
    lambdaU = 10;
    lambdaM = 10;

   
    %% Buffers state at the beginning
    BufferUS = zeros(1,N);
    BufferUAP = 0;
    BufferMAP = 0;

    %% Backoff state r=nan means no packet at the buffer
    r = NaN(1,N);           % Backoff counter
    CW = NaN(1,N);          % Contension window

    DUS_start = NaN(1,N);   % Start time of a unicast packets (for USs), this one used to compute the delay
    
    rM=nan;         % backoff counter of multicast packets
    VL=nan;         % the Value of L: at the beginning is set to nan 
    
    rU=nan;         % The backoff counter of unicast frames of the AP
    CWUAP=nan;      % Contension window of unicast frames of the AP

    DUAP_start=nan; % Start time of a unicast packets (for AP), this one used to compute the delay
    
    CWM = CW_min;   % Contension window of multicast frames, !This one is not doubled when there is a collision!
    %% Initialisations
    US_succ_packet = 0;
    US_lost_packet = 0;
    US_attempt = 0;

    UAP_succ_packet = 0;
    UAP_lost_packet = 0;
    UAP_attempt = 0;

    M_succ_packet = 0;
    M_lost_packet = 0;
    M_attempt = 0;
    
    k=1;
    s=1;
    ss=1;
    
    z=1;

     
    kk1=1; % Time index for THUS for steady state show
    kk2=1; % Time index for LossUS for steady state show

    kk3=1; % Time index for THM
    kk4=1; % Time index for lossM 


%     n = [1:1:100];
%     fprintf('Progress:\n');
%     fprintf(['\n' repmat('.',1,length(n)) '\n\n']);
%     tt = 1;

    current_time = sigma;
%%%%%% Start simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while current_time < T

        %% The probability that a packet is available in the queue
        p_us = 1-exp(-lambdaUS);
        p_uap = 1-exp(-lambdaU);
        p_map = 1-exp(-lambdaM);
        % lambdaUS, lambdaU, lambdaM are the arrival rate of packets to each queue

        %% Buffer check: empty or not?
        % For N USs
        for i=1:N
            if (BufferUS(i) == 0)
                BufferUS(i) = binornd(1,p_us);
                if (BufferUS(i) == 1)
                    CW(i) = CW_min;
                    r(i) = randi([0,CW(i)-1]);
                    DUS_start(i) = current_time;
                else
                    CW(i) = nan;
                    r(i) = nan;
                    DUS_start(i) = nan;
                end
            end
        end
        % For Unicast queue of the AP
        if (BufferUAP == 0)
            BufferUAP = binornd(1,p_uap);
            if (BufferUAP == 1)
                CWUAP = CW_min;
                rU = randi([0,CWUAP-1]);
                DUAP_start = current_time;
            else
                CWUAP = nan;
                rU = nan;
                DUAP_start = nan;
            end 
        end
        % For Multicast queue of the AP
        if (BufferMAP == 0)
            BufferMAP = binornd(1,p_map);
            if (BufferMAP == 1)
                rM = randi([0,CWM-1]);
                DM_start = current_time;
                if (rM == 0)
                    VL = L;
                elseif (rM > 0)
                    VL = nan; 
                end
            else
                rM = nan;
                VL = nan;
                DM_start = nan;
            end
        end

        
        %% check how many USs are ready to transmit
        II = find(r==0); % II is the indexes of all USs that have backoff counter = 0

        Event_flag=0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Find the occurred event

        %% 1-Event UX (US): here, we have US success transmission that is not followed by any multicast packet
        if (length(II)==1) && (rU~=0) && (rM~=0) && (isnan(VL)) && (Event_flag==0)
            
            current_time = current_time + TUX;
            US_succ_packet = US_succ_packet + 1; 
            US_attempt = US_attempt + 1;

            THUS_live(kk1) = US_succ_packet/US_attempt;
            TimeTHUS(kk1) = current_time - TUX;
            kk1 = kk1 + 1;

            BufferUS(II) = BufferUS(II) - 1;
            CW(II)=nan;
            r(II)=nan;

            for i=1:N
                if (i~=II) && (r(i)>0) % all the USs that have a non empty buffer decrease their backoff value by 1
                    r(i) = r(i) - 1;
                end
            end

            if (rU>0)
                rU = rU - 1;
            end

            if (rM>0)
                rM = rM - 1;
                if (rM==0)
                    VL=L;
                end
            end
            
            Delay_US(s) = current_time - DUS_start(II);
                TimeDUS(s) = current_time - TUX;
            s=s+1;
            DUS_start(II) = nan;
               



            Event_flag=1;
        end

        %% 2-Event UX (UAP): here, we have unicast successful transmission by the AP that is not followed by any multicast packet
        if (isempty(II)) && (rU==0) && (isnan(VL)) && (Event_flag==0)

            current_time = current_time + TUX;
            UAP_succ_packet = UAP_succ_packet + 1;
            UAP_attempt = UAP_attempt + 1;

            BufferUAP = BufferUAP - 1;
            CWUAP = nan;
            rU = nan;

            for i=1:N
                if (r(i)>0) % all the USs that have a non empty buffer decrease their backoff value by 1
                    r(i) = r(i) - 1;
                end
            end

            if (rM>0)
                rM = rM - 1;
                if (rM==0)
                    VL = L;
                end
            end
            
            Delay_UAP(ss) = current_time - DUAP_start;
            ss=ss+1;
            clear DUAP_start;

            Event_flag=1;
        end
        %% 3-Event CX  % !here we have a collision between USs!, and no multicast is following that collision
        if (length(II)>1) && (rU~=0) && (rM~=0) && (isnan(VL)) && (Event_flag==0)

            current_time = current_time + TCX;
            US_attempt = US_attempt + length(II);

            for i=1:N
                if (r(i)==0)
                    if (CW(i)==CW_max)
                        US_lost_packet = US_lost_packet + 1;
                            LossUS_live(kk2) = US_lost_packet/US_attempt;
                            Time_Loss_US(kk2) = current_time - TUX;
                            kk2 = kk2 + 1;
                        BufferUS(i) = BufferUS(i) - 1;
                        CW(i) = nan;
                        r(i) = nan;
                        DUS_start(i) = nan;
                    elseif (CW(i)<CW_max)
                        CW(i) = CW(i)*2;
                        r(i) = randi([0,CW(i)-1]);
                    end

                elseif (r(i)>0)
                    r(i) = r(i)-1;
                end
            end

            if (rU>0)
                rU = rU - 1;
            end

            if (rM>0)
                rM = rM - 1;
                if (rM==0)
                    VL=L;
                end
            end

            Event_flag=1;
        end
        %% 4-Event CX  % !here we have a collision between UAP and USs, and no multicast is following that collision
        if (length(II)>=1) && (isnan(VL)) && (rU==0) && (Event_flag==0)

            current_time = current_time + TCX;
            US_attempt = US_attempt + length(II);
            UAP_attempt = UAP_attempt + 1;

            for i=1:N
                if (r(i)==0)
                    if (CW(i)==CW_max)
                        US_lost_packet = US_lost_packet + 1;
                            LossUS_live(kk2) = US_lost_packet/US_attempt;
                            Time_Loss_US(kk2) = current_time - TUX;
                            kk2 = kk2 + 1;
                        BufferUS(i) = BufferUS(i) - 1;
                        CW(i) = nan;
                        r(i) = nan;
                        DUS_start(i) = nan;

                    elseif (CW(i) < CW_max)
                        CW(i) = CW(i)*2;
                        r(i) = randi([0,CW(i)-1]);
                    end

                elseif (r(i)>0)
                    r(i) = r(i) - 1;
                end
            end

            if (CWUAP==CW_max)
                UAP_lost_packet = UAP_lost_packet + 1;
                BufferUAP = BufferUAP - 1;
                CWUAP = nan;
                rU = nan;
                DUAP_start = nan;

            elseif (CWUAP < CW_max)
                CWUAP = 2*CWUAP;
                rU = randi([0,CWUAP-1]);
            end

            if (rM>0)
                rM = rM - 1;
                if (rM==0)
                    VL=L;
                end
            end


            Event_flag=1;
        end

        %% 5-Event collision between unicast (USs) and multicast packets (AP)
        if (length(II)>=1) && (rU~=0) && (VL==0) && (Event_flag==0)

            current_time = current_time + TCX;
            US_attempt = US_attempt + length(II);
            M_attempt = M_attempt + 1;
            M_lost_packet = M_lost_packet + 1;
                LossM_live(kk4) = M_lost_packet/M_attempt;
                Time_LossM(kk4) = current_time - TCX;
                kk4 = kk4 + 1;

            clear DM_start;

            for i=1:N
                if (r(i)==0)
                    if (CW(i)==CW_max)
                        US_lost_packet = US_lost_packet + 1;
                            LossUS_live(kk2) = US_lost_packet/US_attempt;
                            Time_Loss_US(kk2) = current_time - TUX;
                            kk2 = kk2 + 1;
                        BufferUS(i) = BufferUS(i) - 1;
                        CW(i) = nan;
                        r(i) = nan;
                        DUS_start(i) = nan;

                    elseif (CW(i) < CW_max)
                        CW(i) = CW(i)*2;
                        r(i) = randi([0,CW(i)-1]);
                    end
                elseif (r(i)>0)
                    r(i) = r(i) - 1;
                end
            end

            if (rU>0)
                rU = rU - 1;
            end

            BufferMAP = BufferMAP - 1;
            rM = nan;
            VL = nan;
              
            Event_flag=1;
        end

        %% 6-Event UM             unicast packet (of US) followed by a Multicast packet
        if (length(II)==1) && (rU~=0) && (VL>0) && (Event_flag==0)

            current_time = current_time + TUM;
            US_succ_packet = US_succ_packet + 1;
            M_succ_packet = M_succ_packet + 1;
            M_attempt = M_attempt + 1;
                THM_live(kk3) = M_succ_packet/M_attempt;
                Time_THM(kk3) = current_time - (DATA+DIFS+delta);
                kk3 = kk3 + 1;

            US_attempt = US_attempt + 1;

            THUS_live(kk1) = US_succ_packet/US_attempt;
            TimeTHUS(kk1) = current_time - TUM;
            kk1 = kk1 + 1;


            for i=1:N
                if (i~=II) && (r(i)>0)
                    r(i) = r(i) - 1;
                end
            end

            if (rU>0)
                rU = rU - 1;
            end

            BufferUS(II) = BufferUS(II) - 1;
            CW(II) = nan;
            r(II) = nan;

            BufferMAP = BufferMAP - 1;
            rM = nan;
            VL = nan;
            
            Delay_M(k) = current_time - DM_start;
            k=k+1;
            clear DM_start     
            
            Delay_US(s) = current_time - (DATA+DIFS+delta) - DUS_start(II);
                TimeDUS(s) = current_time - TUM;
            s=s+1;
            DUS_start(II) = nan;
                
            
            arrival_time(z) = current_time;
            z=z+1;

            Event_flag=1;
        end

        %% 7-Event UM       unicast packet (of AP) followed by a Multicast packet
        if (isempty(II)) && (rU==0) && (VL>=0) && (Event_flag==0)

            current_time = current_time + TUM;
            UAP_succ_packet = UAP_succ_packet + 1;
            M_succ_packet = M_succ_packet + 1;
            M_attempt = M_attempt + 1;
                THM_live(kk3) = M_succ_packet/M_attempt;
                Time_THM(kk3) = current_time - (DATA+DIFS+delta);
                kk3 = kk3 + 1;
            UAP_attempt = UAP_attempt + 1;

            for i=1:N
                if (r(i)>0)
                    r(i) = r(i) - 1;
                end
            end

            BufferUAP = BufferUAP - 1;
            rU = nan;
            CWUAP = nan;

            BufferMAP = BufferMAP - 1;
            rM = nan;
            VL = nan;
            
            Delay_M(k) = current_time - DM_start;
            k=k+1;
            clear DM_start

            Delay_UAP(ss) = current_time - (DATA+DIFS+delta) - DUAP_start;
            ss=ss+1;
            clear DUAP_start;           
            
            arrival_time(z) = current_time;
            z=z+1;

            Event_flag=1;
        end

        %% 8-Event CM  % here we have collision between USs followed by a multicast packet
        if (length(II)>1) && (rU~=0) && (VL>0) && (Event_flag==0)

            current_time = current_time + TCM;
            M_succ_packet = M_succ_packet + 1;
            M_attempt = M_attempt + 1;
                THM_live(kk3) = M_succ_packet/M_attempt;
                Time_THM(kk3) = current_time - (DATA+DIFS+delta);
                kk3 = kk3 + 1;
            US_attempt = US_attempt + length(II);

            for i=1:N
                if (r(i)==0)
                    if (CW(i)==CW_max)
                        US_lost_packet = US_lost_packet + 1;
                            LossUS_live(kk2) = US_lost_packet/US_attempt;
                            Time_Loss_US(kk2) = current_time - TUX;
                            kk2 = kk2 + 1;
                        BufferUS(i) = BufferUS(i) - 1;
                        CW(i) = nan;
                        r(i) = nan;
                        DUS_start(i) = nan;

                    elseif (CW(i) < CW_max)
                        CW(i) = CW(i)*2;
                        r(i) = randi([0,CW(i)-1]);
                    end

                elseif (r(i)>0)
                    r(i) = r(i) - 1;
                end
            end

            if (rU>0)
                rU = rU - 1;
            end

            BufferMAP = BufferMAP - 1;
            rM = nan;
            VL = nan;

            Delay_M(k) = current_time - DM_start;
            k=k+1;
            clear DM_start;

            arrival_time(z) = current_time;
            z=z+1;

            Event_flag=1;
        end

        %% 9-Event CM  % here we have a collision between USs and UAP, followed by a multicast packet
        if (length(II)>=1) && (rU==0) && (VL>=0) && (Event_flag==0)

            current_time = current_time + TCM;
            M_succ_packet = M_succ_packet + 1;
            M_attempt = M_attempt + 1;
                THM_live(kk3) = M_succ_packet/M_attempt;
                Time_THM(kk3) = current_time - (DATA+DIFS+delta);
                kk3 = kk3 + 1;
            US_attempt = US_attempt + length(II);
            UAP_attempt = UAP_attempt + 1;

            for i=1:N
                if (r(i)==0)
                    if (CW(i)==CW_max)
                        US_lost_packet = US_lost_packet + 1;
                            LossUS_live(kk2) = US_lost_packet/US_attempt;
                            Time_Loss_US(kk2) = current_time - TUX;
                            kk2 = kk2 + 1;
                        BufferUS(i) = BufferUS(i) - 1;
                        CW(i) = nan;
                        r(i) = nan;
                        DUS_start(i) = nan;

                    elseif (CW(i) < CW_max)
                        CW(i) = CW(i)*2;
                        r(i) = randi([0,CW(i)-1]);
                    end

                elseif (r(i)>0)
                    r(i) = r(i) - 1;
                end
            end

            BufferMAP = BufferMAP - 1;
            rM = nan;
            VL = nan;
            
            Delay_M(k) = current_time - DM_start;
            k=k+1;
            clear DM_start;


            if (CWUAP==CW_max)
                UAP_lost_packet = UAP_lost_packet + 1;
                BufferUAP = BufferUAP - 1;
                CWUAP = nan;
                rU = nan;
                DUAP_start = nan;

            elseif (CWUAP < CW_max)
                CWUAP=2*CWUAP;
                rU = randi([0,CWUAP-1]);
            end

            arrival_time(z) = current_time;
            z=z+1;

            Event_flag=1;
        end

        %% 10-Event XM: here when the AP transmit the multicast packet after the expiration of L
        if (isempty(II)) && (rU~=0) && (VL==0) && (Event_flag==0)

            current_time = current_time + TXM;
            M_succ_packet = M_succ_packet + 1;
            M_attempt = M_attempt + 1;
                THM_live(kk3) = M_succ_packet/M_attempt;
                Time_THM(kk3) = current_time - TXM;
                kk3 = kk3 + 1;

            for i=1:N
                if (r(i)>0)
                    r(i) = r(i) - 1;
                end
            end

            if (rU>0)
                rU = rU - 1;
            end

            BufferMAP = BufferMAP - 1;
            rM = nan;
            VL = nan;
            
            Delay_M(k) = current_time - DM_start;
            k=k+1;
            clear DM_start;

            arrival_time(z) = current_time;
            z=z+1;

            Event_flag=1;
        end

        %% 11-Event idle: when no one transmit, idle slot
        if (isempty(II)) && (rU~=0) && (VL~=0) && (Event_flag==0)

            current_time = current_time + sigma;

            for i=1:N
                if (r(i)>0)
                    r(i) = r(i) - 1;
                end
            end

            if (rU>0)
                rU = rU - 1;
            end
            
            if (VL>0)
                VL = VL - 1;
            elseif (rM>0)
                rM = rM - 1;
                if (rM==0)
                    VL=L;
                end
            end

            Event_flag=1;
        end



        %% If no event happened: This one to check if everything is OK
        if (Event_flag==0)
            disp('Stoped because No Event is found');
            r
            rU
            rM
            II
            VL
            return;
        end

        % show the progress of the simulation
%         if (current_time>= (tt/100)*T)
%             fprintf('\b|\n');
%             tt = tt + 1;
%         end
    end % end while
    
    
    
%% The results
    % Throughput
        THUAP = UAP_succ_packet / UAP_attempt;
        THM = M_succ_packet / M_attempt;
        THUS = US_succ_packet / US_attempt;

    % Loss rate
        PUAPloss = UAP_lost_packet / UAP_attempt;
        PMloss = M_lost_packet /M_attempt;
        PUSloss = US_lost_packet / US_attempt;
        
    % Delay
        if exist('Delay_M','var') == 1
            DM = mean(Delay_M);
        else
            DM = nan;
        end

        if exist('Delay_US','var') == 1
            DUS = mean(Delay_US);
        else
            DUS = nan;
        end

        if exist('Delay_UAP','var') == 1
            DUAP = mean(Delay_UAP);
        else
            DUAP = nan;
        end

    % jitter of multicast packets    
        La=length(arrival_time);
        inter_arrival = zeros(1,La);
        inter_arrival(1) = arrival_time(1);

        for k=2:La
            inter_arrival(k) = arrival_time(k) - arrival_time(k-1);
        end
      
        % Jitter = The standard deviation of latency
        JtM = sqrt(var(inter_arrival,1));

end  % end function: MCP_NonSat(N,L)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%