% MakeLightOdour.m
% this script makes kontrol paradigms for an experiment where we present
% both light and odor flickering stimuli to a ORN. It makes a variety of
% kontrol paradigm files:
% 
% 1. (log) steps of light. no odor.
% 2. constant odor, log flickering light
% 3. constant odor, flickering light
% 4. constant light, flickering odor
%
% started 19/1/2015 by srinivas.gs
%
% this script assumes that there is a main air at 2L/min controlled by a
% digital switch, and one additional MFC that is rapidly varied to achive
% the necessary flicker. 
% 
% it uses a frozen noise random stream object called rand_state to get a
% reproducible random numbers (for future experiments)
%
% ControlParadigm Order (should match setup, and config. in Kontroller)
%
% 1. (AO)   to LED Dimmer (0-10V)
% 2. (AO)   to MFC200 (0-5V)
% 2. (DO)   to switch controlling main air @ 2L/min


% some global parameters (defaults)
tc_odour = .1; % switching time. a smaller switching time (like 50ms) is too fast for the MFC to follow
tc_light = .01; 
baseline_dilution = .1/100;
main_flow = 2000; %ml/min
MFC_scale = 40; % 1V = 40mL/min
s = [0 .1 .3 .5 1 1.25 1.5 2 3]; % standard deviation of noise, in units of dilution (%)
dt = 1e-4;



%% section 1 light steps

T = 10;
light_on = floor(1/dt);
light_off = floor(2/dt);

% make initial startup
clear ControlParadigm
ControlParadigm(1).Name = 'CleanAir';
ControlParadigm(1).Outputs = zeros(3,1e4);
ControlParadigm(1).Outputs(3,:) = 1;

% make light steps
light_steps=logspace(log10(1.7),1,10);
for i = 1:length(light_steps)
    ControlParadigm(i+1).Name = strcat('Light-',oval(light_steps(i),2),'V');
    ControlParadigm(i+1).Outputs = zeros(3,floor(T/dt));
    ControlParadigm(i+1).Outputs(3,:) = 1;
    ControlParadigm(i+1).Outputs(1,light_on:light_off) = light_steps(i);
end

%% section 1.1 long light steps

T = 60;
light_on = floor(5/dt);
light_off = floor(55/dt);

% make initial startup

ControlParadigm2(1).Name = 'CleanAir';
ControlParadigm2(1).Outputs = zeros(3,1e4);
ControlParadigm2(1).Outputs(3,:) = 1;

% make light steps
light_steps=logspace(log10(1.7),1,10);
for i = 1:length(light_steps)
    ControlParadigm2(i+1).Name = strcat('Light-',oval(light_steps(i),2),'V');
    ControlParadigm2(i+1).Outputs = zeros(3,floor(T/dt));
    ControlParadigm2(i+1).Outputs(3,:) = 1;
    ControlParadigm2(i+1).Outputs(1,light_on:light_off) = light_steps(i);
end

n = length(ControlParadigm2)+1;
ControlParadigm2(n).Name = 'AllOFF';
ControlParadigm2(n).Outputs = zeros(3,floor(1/dt));

ControlParadigm= [ControlParadigm ControlParadigm2];

% save it
n = ('c:\srinivas\LightSteps_Kontroller_paradigm.mat');

save(n,'ControlParadigm')


