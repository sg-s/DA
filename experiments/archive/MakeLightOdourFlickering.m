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
% this specific script makes one odour stimulus, and combines it with
% differnet levels of light activation. 
%
% ControlParadigm Order (should match setup, and config. in Kontroller)
%
% 1. (AO)   to LED Dimmer (0-5V)
% 2. (AO)   to Microcope light. always zero.
% 3. (AO)   to MFC500 (0-5V)
% 4. (DO)   to switch controlling main air @ 2L/min

clear ControlParadigm
ControlParadigm(1).Name = 'CleanAir';
ControlParadigm(1).Outputs = zeros(4,1e4);
ControlParadigm(1).Outputs(4,:) = 1;


dt  =1e-4;
T = 60;

%switching time
tc = .05; % 50ms is too fast for the MFCs to follow

nsteps= T/tc;
% noise = rand(nsteps,1);
% load frozen noise
load('noise.mat')
noise = [noise noise];
load('gauss_noise.mat')
gauss_noise = [gauss_noise gauss_noise];

% define light levels
light_levels = [0 .8 1 1.5 2 4];
light_gauss_mean = 2.5;
light_gauss_min = .5;
light_gauss_max = 5;

% define light oscillation
light_w = 100; % Hz
light_tau = 1e4/light_w;

% define odour levels
odour_levels = [0 0.01 0.02 0.03 0.05 .1 .5 1];

% for MFC500
dil_min  = 0/100;
dil_max  = .75/100;

dil_mean = (dil_min +dil_max)/2;
V = 2000; %mL/min
MFC_Scale = 100; % 1V=100mL/min
min_flow = 5/200; % guesstimate from turn down ratio

% initialise control paradigm
clear ControlParadigm
ControlParadigm.Name = '';
ControlParadigm.Outputs = [];
c = 1;

% make a clean air paradigm
ControlParadigm.Name = 'CleanAir';
ControlParadigm.Outputs = zeros(4,1e4);
ControlParadigm.Outputs(4,:) = 1;
c = 2;


%% make ShortPulses (light)
for i = 1:length(light_levels)
    ControlParadigm(c).Name = strcat('ShortPulses_',oval(light_levels(i),2),'V');
    ControlParadigm(c).Outputs = zeros(4,1e5); % 10 seconds
    ControlParadigm(c).Outputs(1,4e4:5e4) = light_levels(i);
    % turn main air on
    ControlParadigm(c).Outputs(4,:) = 1;
    c=c+1;
end

%% make ShortPulses (odour)
for i = 1:length(odour_levels)
    ControlParadigm(c).Name = strcat('ShortPulsesOdour_',oval(odour_levels(i),2),'V');
    ControlParadigm(c).Outputs = zeros(4,1e5); % 10 seconds
    ControlParadigm(c).Outputs(3,4e4:5e4) = odour_levels(i);
    % turn main air on
    ControlParadigm(c).Outputs(4,:) = 1;
    c=c+1;
end

%% make long Pulses (light)
for i = 1:length(light_levels)
    ControlParadigm(c).Name = strcat('LongPulses_',oval(light_levels(i),2),'V');
    ControlParadigm(c).Outputs = zeros(4,6e5); % 60 seconds
    % ControlParadigm(c).Outputs(1,5e4:55e4) = light_levels(i);
    
    temp =repmat([ones(light_tau/2,1); zeros(light_tau/2,1)],(T*1e4)/light_tau,1);
    temp(temp==1) = light_levels(i);
    ControlParadigm(c).Outputs(1,:) = temp;
    ControlParadigm(c).Outputs(1,end) = 0;
    ControlParadigm(c).Outputs(1,1:5e4) = 0;
    
    % turn main air on
    ControlParadigm(c).Outputs(4,:) = 1;
    c=c+1;
end


%% make the odour flicker + light background
for i = 1:length(light_levels)
    n = strcat('OdourFlicker+',oval(light_levels(i),2),'V Light');
    ControlParadigm(c).Name = n;
    
    ControlParadigm(c).Outputs= zeros(4,T/dt);
    
    % set to a random number every tc
    a = 1;
    z = a + tc/dt;
    for j = 1:nsteps
        this_dil = noise(j)*(dil_max - dil_min) + dil_min;
        
        
        ControlParadigm(c).Outputs(3,a:z) = (this_dil*V/(1 - this_dil))/MFC_Scale;
        
        % this is for direct control of the flicker amplitude
        % ControlParadigm(c).Outputs(2,a:z) = noise(j)*(odour_levels(end)-odour_levels(2));
        
        % increment
        a = z;
        z = min([a + tc/dt T/dt]);
    end
    
    
    m = mean(ControlParadigm(c).Outputs(3,:));
    ControlParadigm(c).Outputs(3,end) = m; % don't leave the ORN with some high or low value
    ControlParadigm(c).Outputs(3,ControlParadigm(c).Outputs(3,:) > 5) = 5; % clip for sanity
    ControlParadigm(c).Outputs(3,ControlParadigm(c).Outputs(3,:) < min_flow) = min_flow; % clip for sanity

    % add the main air
    ControlParadigm(c).Outputs(4,:) = 1;
    
    % add the light
    % ControlParadigm(c).Outputs(1,:) = light_levels(i);
    temp =repmat([ones(light_tau/2,1); zeros(light_tau/2,1)],(T*1e4)/light_tau,1);
    temp(temp==1) = light_levels(i);
    ControlParadigm(c).Outputs(1,:) = temp;
    ControlParadigm(c).Outputs(1,end) = 0;
    
    c = c+1;
    
end

%% make the light flicker

% fix the light flicker gaussian noise to the correct mean
gauss_noise = 1.5*gauss_noise + light_gauss_mean;

for i = 1:length(odour_levels)
    n = strcat('LightFlicker+',oval(odour_levels(i),2),'V Odour');
    ControlParadigm(c).Name = n;
    
    ControlParadigm(c).Outputs= zeros(4,T/dt);
    
    % set to a random number every tc
    a = 1;
    z = a + tc/dt;
    for j = 1:nsteps

        
        ControlParadigm(c).Outputs(1,a:z) = gauss_noise(j);
        
        % increment
        a = z;
        z = min([a + tc/dt T/dt]);
    end
   
    % smooth the transitions
    led = (ControlParadigm(c).Outputs(1,:));
    transitions = find(diff(led));
    smooth_width = 100;
    for k = 1:length(transitions)
        a = transitions(k) - smooth_width;
        z = transitions(k) + smooth_width;
        ss = (led(z) - led(a))/(z-a);
        led(a:z) = led(a):ss:led(z);
    end
    
    
    ControlParadigm(c).Outputs(1,:) = led;
    
    m = mean(ControlParadigm(c).Outputs(1,:));
    ControlParadigm(c).Outputs(1,end) = m; % don't leave the ORN with some high or low value
    ControlParadigm(c).Outputs(1,ControlParadigm(c).Outputs(1,:) > light_gauss_max) = light_gauss_max; % clip for sanity
    ControlParadigm(c).Outputs(1,ControlParadigm(c).Outputs(1,:) < light_gauss_min) = light_gauss_min; % clip for sanity

    % superimpose a 100Hz 50% duty cycle flicker on top of this
%     blank = repmat([zeros(50,1); ones(50,1)],(T/dt)/100,1);
%     ControlParadigm(c).Outputs(1,:) = ControlParadigm(c).Outputs(1,:).*blank';
%     
    % add the main air
    ControlParadigm(c).Outputs(4,:) = 1;
    
    % add the odour
    ControlParadigm(c).Outputs(3,:) = odour_levels(i);
    
    c = c+1;
    
end


%% add an end
ControlParadigm(c).Name  = 'end';
ControlParadigm(c).Outputs = zeros(4,1e4);

%% save
n = ('c:\srinivas\LightOdourFlicker_Kontroller_paradigm.mat');
save(n,'ControlParadigm')


