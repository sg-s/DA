% MakeShiftedGaussians.m
% this script makes kontrol paradigms for an experiment where we present
% both light and odor flickering stimuli to a ORN. it assumes there are two
% MFCs. MFC500 makes the flickering odor, and MFC200 sets the background level.  
% 
% thus, the flickering sitmulus is always made by MFC500, and the mean of
% the signal is either set by increasing the background level (using
% MFC20) or by adding light. 
%
% started 17/2/2015 by srinivas.gs
%
% it uses a frozen noise random stream object called rand_state to get a
% reproducible random numbers (for future experiments)
%
% ControlParadigm Order (should match setup, and config. in Kontroller)
%
% 1. (AO)   to LED Dimmer (0-10V)
% 2. (AO)   to MFC200 (0-5V)
% 3. (AO)   to MFC500 (0-5V)
% 2. (DO)   to switch controlling main air @ 2L/min


% some global parameters (defaults)
tc_odour = .1; % switching time. a smaller switching time (like 50ms) is too fast for the MFC to follow
tc_light = .01; 
main_flow = 2000; %ml/min
MFC200_scale = 40; % 1V = 40mL/min
MFC500_scale = 100; % 1V = 100mL/min
dt = 1e-4;
T = 120;
baseline_dilution = [0 4.5 7 8 9]/100; % in % dilution, excluding flickering odor and flow

% flicker parmaeters
min_dil = 0; 
max_dil = 7/100; % in % dilution, excluding baseline odor and flow

%% make some startup paradigms
clear ControlParadigm
ControlParadigm(1).Name = 'CleanAir';
ControlParadigm(1).Outputs = zeros(4,1e4);
ControlParadigm(1).Outputs(4,:) = 1;

baseline_V = ((main_flow*baseline_dilution(2))/(1-baseline_dilution(2)))/MFC200_scale;
ControlParadigm(2).Name = 'Baseline';
ControlParadigm(2).Outputs = zeros(4,1e4);
ControlParadigm(2).Outputs(4,:) = 1;
ControlParadigm(2).Outputs(2,:) = baseline_V;

%% make a flicker paradigm

for i = 1:length(baseline_dilution)

    add_here = length(ControlParadigm)+1;

    ControlParadigm(add_here).Name = strcat('Flicker+',oval(100*baseline_dilution(i),1),'%');
    ControlParadigm(add_here).Outputs = zeros(4,T/dt);
    ControlParadigm(add_here).Outputs(4,:) = 1;

    stream = RandStream.getGlobalStream;
    load rand_state
    stream.State = savedState;
    noise = randi(100,(T/tc_odour),1);
    levels = min_dil:(((max_dil-min_dil)/100)):max_dil;

           
    baseline_V = ((main_flow*baseline_dilution(i))/(1-baseline_dilution(i)))/MFC200_scale;
    MFC200_flow = baseline_V*MFC200_scale;
    total_flow = main_flow + MFC200_flow;
    
    a = 1;
    z = a + floor(tc_odour/dt);
    nsteps = T/tc_odour;
    for j = 1:nsteps
        this_dil = levels(noise(j));
        this_V = ((total_flow*this_dil)/(1-this_dil))/MFC500_scale;
        ControlParadigm(add_here).Outputs(3,a:z) = this_V;

        % increment
        a = z;
        z = min([T/dt a+tc_odour/dt]);
    end
    ControlParadigm(add_here).Outputs(3,end) = mean(ControlParadigm(add_here).Outputs(3,:)); % don't leave the ORN with some high or low value
    ControlParadigm(add_here).Outputs(3,ControlParadigm(add_here).Outputs(3,:) > 5) = 5; % clip for sanity
    
    % also led
    ControlParadigm(add_here).Outputs(1,:) = 10*ControlParadigm(add_here).Outputs(3,:);
    ControlParadigm(add_here).Outputs(1,ControlParadigm(add_here).Outputs(1,:) > 10) = 10; % clip for sanity
    ControlParadigm(add_here).Outputs(1,end)=  0; 

    ControlParadigm(add_here).Outputs(2,:) = baseline_V;


end



%% make a close all paradigm
n = length(ControlParadigm)+1;
ControlParadigm(n).Name = 'AllOFF';
ControlParadigm(n).Outputs = zeros(4,floor(1/dt));


%% save
n = ('c:\srinivas\GaussianFlicker_Kontroller_paradigm.mat');
save(n,'ControlParadigm')

