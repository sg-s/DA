% makes control paradigms so that the effective flow through the odor in the
% main air stream is uniformly distrubted, with mean $dil_mean and minimum
% $dil_mean - $dil_rande and maximum $dil_mean + $dil_range
% 
% this is the script used to make the control paradigms in the first mean
% shifted experiment. it has been slighlty modified: to a slower correlaton
% time (100ms), and now there is no step on to a flicker -- it always
% flickers. also, because the configutation has changed, some parts have
% been modified to talk to the MFCs correctly. 
%
% ControlParadigm Order (should match setup, and config. in Kontroller)
%
% 1. (AO)   to LED Dimmer (0-10V)
% 2. (AO)   to MFC200 (0-5V)
% 3. (AO)   to MFC500 (0-5V)
% 2. (DO)   to switch controlling main air @ 2L/min

clear ControlParadigm
ControlParadigm(1).Name = 'CleanAir';
ControlParadigm(1).Outputs = zeros(4,1e4);
ControlParadigm(1).Outputs(4,:) = 1;


dt  =1e-4;
T = 60;

%switching time
tc = .1; % 50ms is too fast for the MFCs to follow

nsteps= T/tc;
% noise = rand(nsteps,1);
% load frozen noise
load('noise.mat')



% for MFC500
dil_min  = [0 2.5 4  9]/100;
dil_max  = [6 8.5 10 11]/100;

dil_mean = (dil_min +dil_max)/2;
V = 2000; %mL/min
MFC_Scale = 100; % 1V=100mL/min

for i = 1:length(dil_mean)
    
    n = strcat('MFC500-',mat2str(100*dil_mean(i)),'%');
    ControlParadigm2(i).Name = n;
    
    ControlParadigm2(i).Outputs= zeros(4,T/dt);
    
    % set to a random number every tc
 
    a = 1;
    z = a + tc/dt;
    for j = 1:nsteps
        this_dil = noise(j)*(dil_max(i) - dil_min(i)) + dil_min(i);
        
        ControlParadigm2(i).Outputs(3,a:z) = (this_dil*V/(1 - this_dil))/MFC_Scale;
        
        
        % increment
        a = z;
        z = min([a + tc/dt T/dt]);
    end
    
    ControlParadigm2(i).Outputs(3,end) = mean(ControlParadigm2(i).Outputs(3,:)); % don't leave the ORN with some high or low value
    ControlParadigm2(i).Outputs(3,ControlParadigm2(i).Outputs(3,:) > 5) = 5; % clip for sanity
    ControlParadigm2(i).Outputs(3,ControlParadigm2(i).Outputs(3,:) < 0) = 0; % clip for sanity

    % add the main air
    ControlParadigm2(i).Outputs(4,:) = 1;
    
    
end

ControlParadigm = [ControlParadigm ControlParadigm2];

l = length(ControlParadigm)+1;

ControlParadigm(l).Name = 'end';
ControlParadigm(l).Outputs = zeros(4,1e4);

%% save
n = ('c:\srinivas\GaussianFlicker_Kontroller_paradigm.mat');
save(n,'ControlParadigm')

