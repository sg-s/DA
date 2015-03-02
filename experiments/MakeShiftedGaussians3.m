% makes control paradigms based on numerical optimistion in MakingMeanShiftedGaussians.m
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



dil_mean = (dil_min +dil_max)/2;
V = 2000; %mL/min
MFC_Scale = 100; % 1V=100mL/min


p.   mu1= 0.2137;
p.sigma1= 0.0362;
p.   mu2= 0;
p.sigma2= 0.4766;
p.  xmin= 0.0133;
p.  xmax= 0.0667;

i = 2;
    
n = strcat('MFC500-.5V');
ControlParadigm2(i).Name = n;

ControlParadigm2(i).Outputs= zeros(4,T/dt);

[~,ControlParadigm2(i).Outputs(3,:)] = BestDistribution([],p);
    
ControlParadigm2(i).Outputs(3,end) = mean(ControlParadigm2(i).Outputs(3,:)); % don't leave the ORN with some high or low value
ControlParadigm2(i).Outputs(3,ControlParadigm2(i).Outputs(3,:) > 5) = 5; % clip for sanity
ControlParadigm2(i).Outputs(3,ControlParadigm2(i).Outputs(3,:) < 0) = 0; % clip for sanity

% add the main air
ControlParadigm2(i).Outputs(4,:) = 1;


ControlParadigm = [ControlParadigm ControlParadigm2];

l = length(ControlParadigm)+1;

ControlParadigm(l).Name = 'end';
ControlParadigm(l).Outputs = zeros(4,1e4);

%% save
n = ('c:\srinivas\ShapedGaussianFlicker_Kontroller_paradigm.mat');
save(n,'ControlParadigm')

