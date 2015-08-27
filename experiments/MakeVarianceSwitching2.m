% MakeVarianceSwitching2.m
% similar to MakeVarianceSwitching, but uses 2MFCs and 2 valves to achieve the switch. 
% 
% makes control paradigms for Kontroller where the variance changes between 2 values, while the mean is kept constant, and this cycles over and over again. 
% 
% this script assumes that there is a main air at 2L/min controlled by a
% digital switch, and one additional MFC that is rapidly varied to achieve
% the necessary flicker. 
% 
%
% ControlParadigm Order (should match setup, and config. in Kontroller)
%
% 1. (AO)   to MFC1
% 2. (AO)   to MFC2  
% 3. (DO)   to switch controlling main air @ 2L/min
% 4. (DO)   to valve for MFC1
% 5. (DO)	to valve for MFC2
%
% This script also assumes that the MFC is running the PD algorithm, with parameters P = 2500 and D = 10000. 
% 
% created by Srinivas Gorur-Shandilya at 1:36 , 07 August 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

% some timing parameters
T = 100; % length of trial
dt = 1e-4; % sampling rate
t_switch  = 5; % time of switch, so we present each stimulus for 5 seconds
tau = 2e-2; % 20 ms switching time, this is very fast, and the MFC will not be reproducible. 

n = floor(T/t_switch); % number of epochs
m = floor(t_switch/tau); % individual setpoints / epoch
o = floor(tau/dt); % this is how long a setpoint is held

% dilution parameters
low_var = 2; % percent
high_var = 3; % percent dilution
mean_dil = 3;

MainFlow = 2000; % mL/min
OdourFlow = 500; % mL/min
MFC_Scale = 100; % 1V= 100mL/min

% set rand stream
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',1984)); 

% make a bunch of replicates
nrep  =5;

for i = 1:nrep
	ControlParadigm(i).Outputs  = ones(5,floor(T/dt));
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

	ControlParadigm(i).Outputs(1,:) = f_low;

	% now we do the same for the high setpoints
	r = rand(m*n,1);
	r = r - .5;
	r = r*high_var*2 + mean_dil; % in percent of total flow

	% convert into actual flow 
	f_high = r*2000./(100-r);

	% convert into MFC control voltage
	f_high = f_high/MFC_Scale;

	f_high = (repmat(f_high,1,o)');
	f_high = f_high(:);

	ControlParadigm(i).Outputs(2,:) = f_high;

	% turn the valve on only in blocks for each MFC
	valve = repmat([0 1],o*m,1);
	valve = valve(:);
	valve = repmat(valve,length(f_low)/length(valve),1);

	ControlParadigm(i).Outputs(3,:) = valve;
	ControlParadigm(i).Outputs(4,:) = abs(1-valve);


end


% bundle the script that made this ControlParadigm with the controlparadigm
genScript = fileread(strcat(mfilename,'.m'));

save('VarianceSwitching_2MFC_Kontroller_paradigm.mat','ControlParadigm','genScript');




