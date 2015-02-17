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

n = length(ControlParadigm)+1;
ControlParadigm(n).Name = 'AllOFF';
ControlParadigm(n).Outputs = zeros(3,floor(1/dt));

% save it
n = ('c:\srinivas\LightSteps_Kontroller_paradigm.mat');

save(n,'ControlParadigm')
%% section 1.1 long light steps

T = 60;
light_on = floor(5/dt);
light_off = floor(55/dt);

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

n = length(ControlParadigm)+1;
ControlParadigm(n).Name = 'AllOFF';
ControlParadigm(n).Outputs = zeros(3,floor(1/dt));

% save it
n = ('c:\srinivas\LongLightSteps_Kontroller_paradigm.mat');

save(n,'ControlParadigm')

%% section 2 constant odor, log flickering light
T = 60;

% set up frozen noise
stream = RandStream.getGlobalStream;
load rand_state
stream.State = savedState;
noise = randi(100,(T/tc_light),1);

levels = logspace(log10(1.7),1,100);
baseline_V = ((main_flow*baseline_dilution)/(1-baseline_dilution))/MFC_scale;

% make initial startup
clear ControlParadigm
ControlParadigm(1).Name = 'CleanAir';
ControlParadigm(1).Outputs = zeros(3,1e4);
ControlParadigm(1).Outputs(3,:) = 1;

ControlParadigm(2).Name = 'LogFlicker';
ControlParadigm(2).Outputs = zeros(3,T/dt);
ControlParadigm(2).Outputs(3,:) = 1;
ControlParadigm(2).Outputs(2,:) = baseline_V;
a = 1;
z = a + floor(tc_light/dt);
nsteps = T/tc_light;
for j = 1:nsteps
    ControlParadigm(2).Outputs(1,a:z) = levels(noise(j));

    % increment
    a = z;
    z = a + tc_light/dt;
end
ControlParadigm(2).Outputs(2,end) = 0;

n = length(ControlParadigm)+1;
ControlParadigm(n).Name = 'AllOFF';
ControlParadigm(n).Outputs = zeros(3,floor(1/dt));

n = ('c:\srinivas\LightLogFlicker_Kontroller_paradigm.mat');
save(n,'ControlParadigm')


return

% make the random variations
nsteps= T/tc;

for i = 1:length(s)
   
    ControlParadigm2(i).Name = strcat('s=',mat2str(s(i)));
    ControlParadigm2(i).Outputs= zeros(2,T/dt);
    
    a = 1;
    z = a + tc/dt;
    for j = 1:nsteps
        this_dil = exp(s(i)*noise(j)) -1; % the - 1 is a correction factor because exp(0) = 1
        ControlParadigm2(i).Outputs(1,a:z) = baseline_V+ (this_dil*main_flow/(100 - this_dil))/MFC_scale;
        
        % increment
        a = z;
        z = a + tc/dt;
    end
    
    % clip for sanity
    ControlParadigm2(i).Outputs(1,ControlParadigm2(i).Outputs(1,:) < baseline_V) = baseline_V;
    ControlParadigm2(i).Outputs(1,ControlParadigm2(i).Outputs(1,:) > 5) = 5;

    % add the main air
    ControlParadigm2(i).Outputs(2,:) = 1;
    
    ControlParadigm2(i).Outputs(1,end) = baseline_V;
    
end

ControlParadigm = [ControlParadigm ControlParadigm2];
clear ControlParadigm2

l = length(ControlParadigm) + 1;
ControlParadigm(l).Name = 'end';
ControlParadigm(l).Outputs = zeros(2,1e4);

n = strcat('Large_Variance_Flicker_bl_',mat2str(baseline_dilution*100));
n = strcat(n,'_Kontroller_paradigm.mat');

save(n,'ControlParadigm')
