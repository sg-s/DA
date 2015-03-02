% GaussianConstructor.m
% Gaussian Constructor is an interactive script meant to be run on an experimental set up (i.e., Kontroller is needed) and will interactively design control paradigms so that some target is met. Here, we say that the target is the distribution of PID values. 
% 
% created by Srinivas Gorur-Shandilya at 2:44 , 02 March 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

% global parameters
dt  =1e-4;
T = 60;
tc = .1; % 50ms is too fast for the MFCs to follow


clc
disp('Gaussian Constructor starting...')

% check for target definition
if ~exist('target')
	error('Target must be defined in the base workspace')
else
	figure, hold on
	plot(target.px,target.py,'k')

	% Construct a questdlg with three options
	choice = questdlg('Is this the target distribution?', ...
		'Gaussians', ...
		'No','Yes','Yes');
	% Handle response
	switch choice
	    case 'No'
	        error('wrong distribution')
	    case 'Yes'
	    	delete(gcf)
	end
end

mfc_min = input('Enter minimum of MFC control range.');
mfc_max = input('Enter maximum of MFC control range.');

disp('Sampling uniformly from this range, and making a control paradigm...')
ControlParadigm = MakeGaussianFlicker(mfc_min,mfc_max,T,dt,tc);

disp('DONE. Running the experiment...')

h = msgbox('Ready to run experiment?','GaussianConstructor');
uiwait(h);

data = Kontroller('ControlParadigm',ControlParadigm,'RunTheseParadigms',[1 1 1 1 1 2]);
disp('Censoring the first trial...')
data(1).PID(1,:) = [];
data(1).MFC500(1,:) = [];


disp('Finished running experiment. Will now build a model for the delivery system...')


[K_MFC,K_PID,p_hill] = BuildDeliverySystemModel(data,ControlParadigm,2);
disp('DONE')



