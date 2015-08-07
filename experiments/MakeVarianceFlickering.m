% MakeVarianceFlicker.m
% makes control paradigms for Kontroller where the variance changes along the time series
% 
% this script assumes that there is a main air at 2L/min controlled by a
% digital switch, and one additional MFC that is rapidly varied to achive
% the necessary flicker. 
% 
%
% ControlParadigm Order (should match setup, and config. in Kontroller)
%
% 1. (AO)   to MFC 
% 2. (DO)   to switch controlling main air @ 2L/min
% 
% created by Srinivas Gorur-Shandilya at 1:36 , 07 August 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

% some core parameters 
T = 60; % length of trial
dt = 1e-4; % sampling rate
tau1 = 2e-2; % 20ms, switching time of individual steps
tau2 = 1e-1; % 100ms, switching time of variance. 

mean_dil = 4; % 4%
max_dil = 8;
min_dil = 1;
MainFlow = 2000; % mL/min
OdourFlow = 500; % mL/min
MFC_Scale = 100; % 1V= 100mL/min

nsteps = floor(T/tau1);
msteps = floor(T/tau2);
steps_per_epoch = nsteps/msteps;

% set rand stream
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',1984)); 

variances = rand(msteps,1);

ControlParadigm.Outputs  = ones(2,floor(T/dt));
ControlParadigm.Name = 'VarianceFlickering';

a = 1;
z = a + floor(tau1/dt);
for i = 1:msteps
	% choose permitted minimum and maximum dilutions in this epoch
	this_max_dil = mean_dil + variances(i)*(max_dil - mean_dil);
	this_min_dil = mean_dil - variances(i)*(mean_dil - min_dil);

	% choose dilutions for the steps in this epoch
	these_dilutions = rand(steps_per_epoch,1);
	these_dilutions = this_min_dil + these_dilutions*(this_max_dil - this_min_dil);

	for j = 1:steps_per_epoch
		% figure out the control signal for this desired dilution
		this_control_signal = (these_dilutions(j)*MainFlow)/(100-these_dilutions(j));
		ControlParadigm.Outputs(1,a:z) = this_control_signal/MFC_Scale;

		a = z;
		z = z+floor(tau1/dt);

	end
end


% bundle the script that made this ControlParadigm with the controlparadigm
genScript = fileread(strcat(mfilename,'.m'));

save('Variance_Flickering_Kontroller_Paradigm.mat','ControlParadigm','genScript');




