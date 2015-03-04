% Make Gaussian Flicker from optimised distributions
% this m file looks around for mat files named "dist*.mat", loads them
% it expects each mat file to have a structure called p, which contains parmaeters for
% dist_gauss2.m 
% it uses these parmaeters to create a distribution, and sample from it and construct
% a control paradigm for driving a MFC
% 
% created by Srinivas Gorur-Shandilya at 8:51 , 04 March 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

% global parameters
mfc_min = 0;
mfc_max = 5; % this doesn't matter, as we replace the MFC signal
T = 60;
dt = 1e-4;
tc= .1; 

[ControlParadigm] = MakeGaussianFlicker(mfc_min,mfc_max,T,dt,tc);
AllOff = ControlParadigm(end);
ControlParadigm(end) = [];


allfiles = dir('dist*.mat');
for i = 1:length(allfiles)
	load(allfiles(i).name)
	[~,s] = BestDistribution([],p);
	ControlParadigm(end+1) = ControlParadigm(end);
	s(s<5/200) = 5/200;
	ControlParadigm(end).Outputs(2,:) = s;
	ControlParadigm(end).Name = strcat('MFC',oval(mean(s)));
end

ControlParadigm(end+1) = AllOff;
ControlParadigm(1)=[];
save('Optimised_Gaussian_Flicker_Kontroller_Paradigm.mat','ControlParadigm')