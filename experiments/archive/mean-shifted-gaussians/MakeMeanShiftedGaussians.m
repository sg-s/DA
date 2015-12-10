% MakeMeanShiftedGaussians.m
% makes control paradigms for Kontroller
% such that we get many odour distributions with differing means
% but similar variances 
%
% experimental setup:
% 
% Outputs:
% MFC500
% valve for odourised air
% switch for main air
% 
% created by Srinivas Gorur-Shandilya at 11:06 , 30 June 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

% timing parameters
T = 60;
dt = 1e-4;
tc = .05; % switching time
nsteps = T/tc;
flicker_on = 5;
flicker_off = 55;

% MFC parameters
MFC_scale = 100; % 1V = 100mL/min
Main_flow = 2000; %mL/min
MFC_min_V = 5/200; % V

% distribution parameters
dilution_mean = [1.5:.5:7]/100; % in percent
dilution_var = 3/100;

% fix random stream
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);
noise = rand(10000,1);

for i = 1:length(dilution_mean)
	ControlParadigm(i).Name = strcat('Flicker-',oval(100*dilution_mean(i)),'%');
	ControlParadigm(i).Outputs = zeros(3,T/dt);
	ControlParadigm(i).Outputs(3,:) = 1;
    ControlParadigm(i).Outputs(2,flicker_on/dt:flicker_off/dt) = 1;

	% make the flicker
	dil_min = dilution_mean(i) - dilution_var;
	dil_max = dilution_mean(i) + dilution_var;

	a = 1;
	z = (T/dt)/nsteps;
	for j = 1:nsteps
		this_dilution = dil_min + noise(j)*(dil_max - dil_min);
		req_odour_flow = (Main_flow*this_dilution)/(1 - this_dilution); % in mL/min
		ControlParadigm(i).Outputs(1,a:z) = req_odour_flow/MFC_scale;

		a = z;
		z = a + (T/dt)/nsteps;
    end

    % make sure MFC never fully shuts down
    ControlParadigm(i).Outputs(1,ControlParadigm(i).Outputs(1,:) < MFC_min_V) = MFC_min_V;
    
	% don't leave it at a weird place
	ControlParadigm(i).Outputs(1,end) =  mean(ControlParadigm(i).Outputs(1,:));
    
    
    
end


% bundle the script that made this ControlParadigm with the controlparadigm
genScript = fileread(strcat(mfilename,'.m'));

save('MSG2015_Kontroller_Paradigm.mat','ControlParadigm','genScript');










