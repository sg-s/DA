% makeVarianceSwitching_GCamp6.m
% constructs a variance switching paradigm and images Calcium at the same time too.
% 
% makes control paradigms for Kontroller where the variance changes between 2 values, while the mean is kept constant, and this cycles over and over again. 
% 
% this script assumes that there is a main air at 2L/min controlled by a
% digital switch, and one additional MFC that is rapidly varied to achieve
% the necessary flicker. 
% 
% for well-balanced paradigms with very similar means, you need to manually tune this script with a fudge factor. use findFudgeFactor (also included in this folder) to find the best fudge_factor
%
% ControlParadigm Order (should match setup, and config. in Kontroller)
%
% 1. (AO)	to blue LED illumination for GCamp6 excitation 
% 2. (AO)   to MFC1
% 3. (AO)   to MFC2  
% 4. (DO)   to valve for MFC1
% 5. (DO)	to valve for MFC2
% 6. (DO)   to switch controlling main air @ 2L/min
%
% This script also assumes that the MFC is running the PD algorithm, with parameters P = 2500 and D = 10000. You have to set this manually. 
% 
% created by Srinivas Gorur-Shandilya at 1:36 , 07 August 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function ControlParadigm = makeVarianceSwitching(fudge_factor)


% some timing parameters
T = 200; % length of trial
dt = 1e-4; % sampling rate
t_switch  = 5; % time of switch, so we present each stimulus for 5 seconds
tau = 2e-2; % 20 ms switching time, this is very fast, and the MFC will not be reproducible. 

n = floor(T/t_switch); % number of epochs
m = floor(t_switch/tau); % individual setpoints / epoch
o = floor(tau/dt); % this is how long a setpoint is held

% dilution parameters
low_var = 2;  % percent dilution
high_var = 5; % percent dilution 
mean_dil = 3; % percent dilution
if ~nargin
    fudge_factor = .5;
end

MainFlow = 2000; % mL/min
OdourFlow = 500; % mL/min
MFC_Scale = 100; % 1V = 100mL/min

% set rand stream
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',1984)); 

% make a bunch of replicates
nrep = 5;

for i = 1:nrep
	ControlParadigm(i).Outputs  = ones(6,floor(T/dt));
	ControlParadigm(i).Name = ['VarianceSteps_rep_' mat2str(i)];

	% make the low variance MFC control signal
	r = rand(m*n,1);
	r = r - .5;
	r = r*low_var*2 + mean_dil; % in percent of total flow

	% convert into actual flow 
	f_low = r*2000./(100-r);

	% convert into MFC control voltage
	f_low = f_low/MFC_Scale;

	f_low = (repmat(f_low,1,o)');
	f_low = f_low(:);

	ControlParadigm(i).Outputs(2,:) = f_low;

	% now we do the same for the high setpoints
	r = rand(m*n,1);
	r = r - .5;
	r = r*high_var*2 + mean_dil*fudge_factor; % in percent of total flow

	% convert into actual flow 
	f_high = r*2000./(100-r);

	% convert into MFC control voltage
	f_high = f_high/MFC_Scale;

	f_high = (repmat(f_high,1,o)');
	f_high = f_high(:);

	ControlParadigm(i).Outputs(3,:) = f_high;

	% turn the valve on only in blocks for each MFC
	valve = repmat([0 1],o*m,1);
	valve = valve(:);
	valve = repmat(valve,length(f_low)/length(valve),1);

   	ControlParadigm(i).Outputs(4,:) = abs(1-valve);
	ControlParadigm(i).Outputs(5,:) = valve;
    
    % blank first five seconds
    ControlParadigm(i).Outputs(4,1:5/dt)=0;
    ControlParadigm(i).Outputs(5,1:5/dt)=0;
    
    % blank last five seconds
    ControlParadigm(i).Outputs(4,end-5/dt:end)=0;
    ControlParadigm(i).Outputs(5,end-5/dt:end)=0;
    
    % make sure we dont do something stupid
    ControlParadigm(i).Outputs(2,ControlParadigm(i).Outputs(2,:)<0) = 0;
    ControlParadigm(i).Outputs(3,ControlParadigm(i).Outputs(3,:)<0) = 0;
    ControlParadigm(i).Outputs(2,ControlParadigm(i).Outputs(2,:)>5) = 5;
    ControlParadigm(i).Outputs(3,ControlParadigm(i).Outputs(3,:)>5) = 5;

    % turn on the blue LED to image GCamp6 (and to trigger the data)
    ControlParadigm(i).Outputs(1,:) = 1;
    ControlParadigm(i).Outputs(1,end) = 0;
    
end



% bundle the script that made this ControlParadigm with the controlparadigm
genScript = fileread(strcat(mfilename,'.m'));
save('GCamp6_contrast_switching_Kontroller_paradigm.mat','ControlParadigm','genScript');




