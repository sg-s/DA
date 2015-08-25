% MakeVarianceSwitching.m
% makes control paradigms for Kontroller where the variance changes between 2 values, while the mean is kept constant, and this cycles over and over again. 
% 
% this script assumes that there is a main air at 2L/min controlled by a
% digital switch, and one additional MFC that is rapidly varied to achieve
% the necessary flicker. 
% 
%
% ControlParadigm Order (should match setup, and config. in Kontroller)
%
% 1. (AO)   to MFC 
% 2. (DO)   to switch controlling main air @ 2L/min
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
	ControlParadigm(i).Outputs  = ones(2,floor(2*T/dt));
	ControlParadigm(i).Name = ['VarianceSteps_rep_' mat2str(i)];

	% make the low mean setpoints
	r = rand(m*n,1);
	r = r - .5;
	r = r*low_var*2 + mean_dil; % in percent of total flow

	% convert into actual flow 
	f_low = r*2000./(100-r);

	% convert into MFC control voltage
	f_low = f_low/MFC_Scale;

	% now we do the same for the high setpoints
	r = rand(m*n,1);
	r = r - .5;
	r = r*high_var*2 + mean_dil; % in percent of total flow

	% convert into actual flow 
	f_high = r*2000./(100-r);

	% convert into MFC control voltage
	f_high = f_high/MFC_Scale;

	% rearrange the low and high into blocks
	f = [reshape(f_low,length(f_low)/m,m) reshape(f_high,length(f_low)/m,m)];
	f = f'; f  = f(:);

	% expand into full control vector
	f = repmat(f,1,o)';
	f = f(:);

	ControlParadigm(i).Outputs(1,:) = f;
end


% bundle the script that made this ControlParadigm with the controlparadigm
genScript = fileread(strcat(mfilename,'.m'));

save('MeanSwitching.mat','ControlParadigm','genScript');




