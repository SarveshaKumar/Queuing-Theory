% Script/Function Name: MC_Event_D_MM2_Q_V3
% Version: 3.1
% Date: 4/3/18
% Author: Sarvesha Kumar Kombaiah Seetha
% Description: This scritpt uses a MC algorithm to calculate 
% time-dependent waiting times and queue lengths for a simple 
% two-server queueing system. This version doesnt make simplifying assumptions. 
clc 
clf
clear
tic
% Input Module
MIAR = 2;          % Mean interarrival rate (/m)
MIAT = 1/MIAR;       % Mean interarrival time (in m)
MST = 0.667;         % Mean Service time foer each server (in m)
MSR = 1/MST;         % Mean Service Rate (/m)
Nc = 20000;          % Number of customers
Tsys = zeros(1,Nc);    % Create/initialize time in system vector 
Tia = zeros(1,Nc);     % Create/initialize customer interarrival vector 
Ta = zeros(1,Nc);      % Create/initialize customer arrival vector
Tsrv = zeros(1,Nc);    % Create/initialize service time vector 
Tnsrv = zeros(1,Nc);   % Create/initialize enter service time vector
Txsys = zeros(1,Nc);   % Create/initialize exit service time vector 
WTqmc = zeros(1,Nc);   % Create/initialize time in Q vector  
MTqmc = zeros(1,Nc);   % Create/initialize Mean Time in Q vector  
MTsys = zeros(1,Nc);   % Create/initialize Meqan Time in System vector 
% Calculation Module
% Generate 
Timespdat = makedist('exponential', 'mu', MIAT); 
% Generate random arrival times
pdser = makedist('exponential', 'mu', MST); 
% Generate random service times
% Set up using 1st & 2nd customer
Tia(1) = random(pdat);       
% Generate interarrival times for cust 1
Tia(2) = random(pdat);
Ta(1) = Tia(1);             
% Calc arrival time for 1
Ta(2) = Tia(2); Tsys(1) = Tsrv(1);          
% Calc time in system for 1
Tsys(2) = Tsrv(2); Tnsrv(1) = Ta(1);            
% Calc enter service time for 1
Tnsrv(2) = Ta(2);
Tsrv(1) = random(pdser);     
% Generate service time for cust 1
Tsrv(2) = random(pdser);
Txsys(1) = Ta(1) + Tsrv(1);  
% Calc exit system time 
Txsys(2) = Ta(2) + Tsrv(2);
WTqmc(1) = 0;               
% Calc queue waiting time for 1
WTqmc(2) = 0;
MTqmc(1) = WTqmc(1);        
% Initialize Mean time in Queue
MTqmc(2) = WTqmc(2);
% MTsys(1) = Tsys(1);       % Intialize Mean time in System
iInSer = 1; % Used to track customers w. longer service times
% Calc times for remaining customers
for i = 3:1:Nc   
    % Calc arrival time for customer i   
    Tia(i) = random(pdat);       
    Ta(i) = Ta(i-1) + Tia(i);   
    Tsys(i) = Ta(i);    
    % Calc time customer enters service    
    if ((Ta(i) < Txsys(i-1)) && (Ta(i) < Txsys(iInSer)))   
        % Both servers busy when customer i enters.        
        if Txsys(i-1) < Txsys(iInSer)          
            % Customer i-1 leaves service prior to iInSer        
            Tnsrv(i) = Tnsrv(i-1) + Tsrv(i-1);          
            % Customer i enters service when i-1 exits.      
        else  % Customer iInSer leaves before i-1         
            Tnsrv(i) = Tnsrv(iInSer) + Tsrv(iInSer);         
            iInSer = i-1;       
        end
    else      % at least one server open => no queue              
        % => i goes right in       
        Tnsrv(i) = Ta(i);       
        iInSer = i-1;    
    end
    Tsrv(i) = random(pdser);   
    % Calc service time for customer     
    Txsys(i) = Tnsrv(i) + Tsrv(i);  
    % Calc time customer exits system    
    WTqmc(i) = Tnsrv(i)-Ta(i);     
    % Calc time in queue    
    Tsys(i) = Txsys(i) - Ta(i);   
    % Calc time in system    
    MTqmc(i) = ((i-1)* MTqmc(i-1) + WTqmc(i))/i;  
    % Update Mean Time in Queue   
    % MTsys(i) = ((i-1)* MTsys(i-1) + Tsys(i))/i; 
    % Update MT in Systemend
    % Calcuate Queue Length at differnt points in time
    Tend = Txsys(Nc);Nts = floor(Tend);      
    % => plot queue length every hour
    dt = 1;
    Time = zeros(1,Nts);       % Create/initialize System Clock vector.
    Lq = zeros(1,Nts);Time(1) = 0;Lq(1) = 0;                 % Initialize queue
    for j = 2:1:Nts   
        Time(j) = Time(j-1)+dt;    
        Lq(j) = 0;    
        % Find all customers that arrived prior to Time(j), that did not exit   
        % prior to Time j, and have not entered service.    
        for i = 1:1:Nc       
            if (Ta(i) < Time(j) && Txsys(i) > Time(j) && Tnsrv(i) > Time(j))   
                Lq(j) = Lq(j)+1;       
            end
        end
    end
    % Output Module
    rhoa = MST/(2*MIAT)
    MLqa = 2*rhoa^3/((1-rhoa)*(1+rhoa))
    MLqmc = mean(Lq)
MTqa = MIAT * MLqa
MTqmc1 = MTqmc(Nc)
MTqmc2 = mean(WTqmc)
SDmtamc = std(WTqmc)
SEmtamc = SDmtamc/sqrt(Nc)
toc
% plot(MTqmc(t) & Lq(t))
subplot(2,1,1)
plot(Ta,MTqmc)
subplot(2,1,2)
plot(Time,Lq)