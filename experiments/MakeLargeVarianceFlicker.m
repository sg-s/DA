% LargeVarianceFlicker.m
% this script makes kontrol paradigms for an experiment where we present
% randomly flickering odor stimuli to a ORN. the variance of this randomly
% flcikering odor stimuli varies, and we want it to be large so that we can
% look for effects of fast gain changes, if any. 
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
% 1. (AO)   to MFC 
% 2. (DO)   to switch controlling main air @ 2L/min

function [ControlParadigm] = MakeLargeVarianceFlicker(varargin)

% some global parameters (defaults)
tc = .05; % switching time. a smaller switching time (like 50ms) is too fast for the MFC to follow
baseline_dilution = 1/100;
main_flow = 2000; %ml/min
MFC_scale = 40; % 1V = 40mL/min
s = [0 .1 .3 .5 1 1.25 1.5 2 3]; % standard deviation of noise, in units of dilution (%)
T= 60;
dt = 1e-4;

% evaluate optional arguments
for i = 1:nargin
    eval(varargin{i})
end

% set up frozen noise
stream = RandStream.getGlobalStream;
load rand_state
stream.State = savedState;
noise = randn(1,ceil(T/dt));

% make initial startup
clear ControlParadigm
ControlParadigm(1).Name = 'Start';
ControlParadigm(1).Outputs = ones(2,1e4);
baseline_V = ((main_flow*baseline_dilution)/(1-baseline_dilution))/MFC_scale;
ControlParadigm(1).Outputs(1,:) = baseline_V;

% make the random variations
nsteps= T/tc;

for i = 1:length(s)
   
    ControlParadigm2(i).Name = strcat('s=',mat2str(s(i)));
    ControlParadigm2(i).Outputs= zeros(2,T/dt);
    
    a = 1;
    z = a + tc/dt;
    for j = 1:nsteps
        this_dil = exp(s(i)*noise(j));
        
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
    
end

ControlParadigm = [ControlParadigm ControlParadigm2];
clear ControlParadigm2

l = length(ControlParadigm) + 1;
ControlParadigm(l).Name = 'end';
ControlParadigm(l).Outputs = zeros(2,1e4);

save('Large_Variance_Flicker_Kontroller_paradigm.mat','ControlParadigm')
