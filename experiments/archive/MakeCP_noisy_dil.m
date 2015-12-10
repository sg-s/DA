% makes control paradigms so that the effective flow through the odor in the
% main air stream is uniformly distrubted, with mean $dil_mean and minimum
% $dil_mean - $dil_rande and maximum $dil_mean + $dil_range
% 
% this is the script used to make the control paradigms in the first mean
% shifted experiment. 

dil_mean = 0.03;
dil_range = 0.0209; 
V = 2000; %mL/min
MFC_Scale = 40; % 1V=40mL/min

ControlParadigm(1).Name = 'Start';
ControlParadigm(1).Outputs = zeros(5,1e4);
ControlParadigm(1).Outputs(1,:) = 1;
ControlParadigm(1).Outputs(2,:) = 2.5;
ControlParadigm(1).Outputs(5,:) = 1;


dt  =1e-4;
T = 60;
t_on = 5;
t_off = 55;
rand_on = 12;
rand_off = 55;

%switching time
tc = .05; % 50ms

nsteps= (rand_off-rand_on)/tc;
% noise = rand(nsteps,1);
% load frozen noise
load('noise.mat')


% for MFC200
for i = 1:length(dil_mean)
    
    n = strcat('MFC200-',mat2str(100*dil_mean(i)),'%');
    ControlParadigm2(i).Name = n;
    
    ControlParadigm2(i).Outputs= zeros(5,T/dt);
    
    % turn the steady state on
    ControlParadigm2(i).Outputs(2,:) = (dil_mean(i)*V)/MFC_Scale;
    
    % set to a random number every tc
    rmin = dil_mean(i) - dil_range;
    rmax = dil_mean(i) + dil_range;
    
    a = rand_on/dt;
    z = a + tc/dt;
    for j = 1:nsteps
        this_dil = noise(j)*dil_range*2 + rmin;
        
        ControlParadigm2(i).Outputs(2,a:z) = (this_dil*V/(1 - this_dil))/MFC_Scale;
        
        
        % increment
        a = z;
        z = a + tc/dt;
    end

    % add the main air
    ControlParadigm2(i).Outputs(5,:) = 1;

    
    % turn valve for MFC200 on
    ControlParadigm2(i).Outputs(4,t_on/dt:t_off/dt) = 1;
    
    % set the other MFC to some nice value
    ControlParadigm2(i).Outputs(1,:) = 1;
    
    
end

ControlParadigm = [ControlParadigm ControlParadigm2];
clear ControlParadigm2


% for MFC500
dil_mean = [4 5 7]/100; 
dil_range = 0.02; 
V = 2000; %mL/min
MFC_Scale = 100; % 1V=100mL/min

for i = 1:length(dil_mean)
    
    n = strcat('MFC500-',mat2str(100*dil_mean(i)),'%');
    ControlParadigm2(i).Name = n;
    
    ControlParadigm2(i).Outputs= zeros(5,T/dt);
    
    % turn the steady state on
    ControlParadigm2(i).Outputs(1,:) = (dil_mean(i)*V)/MFC_Scale;
    
    % set to a random number every tc
    rmin = dil_mean(i) - dil_range;
    rmax = dil_mean(i) + dil_range;
    
    a = rand_on/dt;
    z = a + tc/dt;
    for j = 1:nsteps
        this_dil = noise(j)*dil_range*2 + rmin;
        
        ControlParadigm2(i).Outputs(1,a:z) = (this_dil*V/(1 - this_dil))/MFC_Scale;
        
        
        % increment
        a = z;
        z = a + tc/dt;
    end

    % add the main air
    ControlParadigm2(i).Outputs(5,:) = 1;
    
    % turn valve for MFC200 on
    ControlParadigm2(i).Outputs(3,t_on/dt:t_off/dt) = 1;
    
    % set the other MFC to some nice value
    ControlParadigm2(i).Outputs(2,:) = 2.5;
    
end

ControlParadigm = [ControlParadigm ControlParadigm2];

l = length(ControlParadigm)+1;

ControlParadigm(l).Name = 'end';
ControlParadigm(l).Outputs = zeros(5,1e4);
