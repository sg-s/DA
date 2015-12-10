% makes control paradigms with long pulses of light with 50% duty cycle and
% diffeent frequencies. the point is to find an optimal ay of stimulating
% ORNs with light for a long period with sustained activity. 
%
% ControlParadigm Order (should match setup, and config. in Kontroller)
%
% 1. (AO)   to LED Dimmer (0-5V)
% 2. (AO)   to MFC500 (0-5V)
% 3. (DO)   to switch controlling main air @ 2L/min

clear ControlParadigm
ControlParadigm(1).Name = 'CleanAir';
ControlParadigm(1).Outputs = zeros(3,1e4);
ControlParadigm(1).Outputs(3,:) = 1;


dt  =1e-4;
T = 30;


% define light levels
light_levels = [1];

% define light switching frequency
w = [20 25 40 50 100 200 500 ]; %Hz
tau = floor(1e4./w);


% initialise control paradigm
clear ControlParadigm
ControlParadigm.Name = '';
ControlParadigm.Outputs = [];
c = 1;

% make a clean air paradigm
ControlParadigm.Name = 'CleanAir';
ControlParadigm.Outputs = zeros(3,1e4);
ControlParadigm.Outputs(3,:) = 1;
c = 2;



%% make long Pulses (light)
for i = 1:length(w)
    ControlParadigm(c).Name = strcat('LightPulses_',oval(w(i),2),'Hz');
    ControlParadigm(c).Outputs = zeros(3,3e5); % 30 seconds
    
    tw = floor(tau(i)/2);
    temp =repmat([ones(tw,1); zeros(tw,1)],(T*1e4)/tau(i),1);
    temp(temp==1) = light_levels;
    ControlParadigm(c).Outputs(1,:) = temp;
    ControlParadigm(c).Outputs(1,1:5e4) = 0;
    ControlParadigm(c).Outputs(1,end) = 0;
    % turn main air on
    ControlParadigm(c).Outputs(3,:) = 1;
    c=c+1;
end

%% add an end
ControlParadigm(c).Name  = 'end';
ControlParadigm(c).Outputs = zeros(3,1e4);

%% save
n = ('c:\srinivas\LightPulses_Kontroller_paradigm.mat');
save(n,'ControlParadigm')

