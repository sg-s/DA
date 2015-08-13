% tunePID
% script to allow you to tune PID settings of MFC
% the target is to get a flicker with a 10ms switching time
% 
% created by Srinivas Gorur-Shandilya at 4:15 , 13 August 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

% some core parameters 
T = 5; % length of trial
dt = 1e-4; % sampling rate
tau1 = 2e-2; % 20ms, switching time of individual steps

mean_dil = 4; % 4%
max_dil = 7;
min_dil = 1;
MainFlow = 2000; % mL/min
OdourFlow = 500; % mL/min
MFC_Scale = 100; % 1V= 100mL/min

nsteps = floor(T/tau1);

% set rand stream
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',1984)); 

ControlParadigm.Outputs  = ones(2,floor(T/dt));
ControlParadigm.Name = 'test';

a = 1;
z = a + floor(tau1/dt);

% choose dilutions for all steps
these_dilutions = rand(nsteps,1);
these_dilutions = min_dil + these_dilutions*(max_dil - min_dil);

for j = 1:nsteps
	% figure out the control signal for this desired dilution
	this_control_signal = (these_dilutions(j)*MainFlow)/(100-these_dilutions(j));
	ControlParadigm.Outputs(1,a:z) = this_control_signal/MFC_Scale;

	a = z;
	z = z+floor(tau1/dt);

end

% make sure it off at the end
ControlParadigm.Outputs(1,end) = 0;

% get the PID values
p = input('enter P value:\n');
i = input('enter I value:\n');
d = input('enter D value:\n');

% run the experiment
data = Kontroller('ControlParadigm',ControlParadigm,'RunTheseParadigms',ones(1,5),'w',1e4);




