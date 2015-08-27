% MakeMeanSwitching2.m
% version of MakeMeanSwitching that uses two MFCs to achieve a fast switch
% 
% makes control paradigms for Kontroller where the mean changes between 2 values, while the variance is kept constant, and this cycles over and over again. 
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
% 3. (DO)   to valve for MFC1
% 4. (DO)   to valve for MFC2
% 5. (DO)	to switch controlling main air @ 2L/min
%
% This script also assumes that the MFC is running the PD algorithm, with parameters P = 2500 and D = 10000. 
% 
% created by Srinivas Gorur-Shandilya at 1:36 , 07 August 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

% some timing parameters
T = 200; % length of trial in seconds
dt = 1e-4; % sampling rate
t_switch  = 5; % time of switch, so we present each stimulus for 5 seconds
tau = 5e-2; % 20 ms switching time, this is very fast, and the MFC will not be reproducible. 

n = floor(T/t_switch); % number of epochs
m = floor(t_switch/tau); % individual setpoints / epoch
o = floor(tau/dt); % this is how long a setpoint is held

% dilution parameters
low_mean = 3; % percent
high_mean = 9; % percent dilution
range = 3; % percent, below and above the mean. 
fudge_factor = 1.5; % scale on the variance for the high setpoint

MainFlow = 2000; % mL/min
OdourFlow = 500; % mL/min
MFC_Scale = 100; % 1V= 100mL/min

% set rand stream
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',1984)); 

% make a bunch of replicates
nrep = 5;

for i = 1:nrep
	ControlParadigm(i).Outputs  = ones(5,floor(T/dt));
	ControlParadigm(i).Name = ['MeanSteps_rep_' mat2str(i)];

	% make the low mean setpoints
	r = rand(m*n,1);
	r = r - .5;
	r = r*range*2 + low_mean; % in percent of total flow

	% convert into actual flow 
	f_low = r*2000./(100-r);

	% convert into MFC control voltage
	f_low = f_low/MFC_Scale;

	f_low = (repmat(f_low,1,o)');
	f_low = f_low(:);
    f_low(f_low<0)= 0;
    f_low(f_low>5) = 5;

	ControlParadigm(i).Outputs(1,:) = f_low;

	% now we do the same for the high setpoints
	r = rand(m*n,1);
	r = r - .5;
	r = r*range*2*fudge_factor + high_mean; % in percent of total flow

	% convert into actual flow 
	f_high = r*2000./(100-r);

	% convert into MFC control voltage
	f_high = f_high/MFC_Scale;

	f_high = (repmat(f_high,1,o)');
	f_high = f_high(:);
    f_high(f_high<0) = 0;
    f_high(f_high>5) = 5;

	ControlParadigm(i).Outputs(2,:) = f_high;


	% turn the valve on only in blocks for each MFC
	valve = repmat([0 1],o*m,1);
	valve = valve(:);
	valve = repmat(valve,length(f_low)/length(valve),1);
    

	ControlParadigm(i).Outputs(4,:) = valve;        
	ControlParadigm(i).Outputs(3,:) = abs(1-valve);
    
    ControlParadigm(i).Outputs(4,1:5e4) = 0;         
	ControlParadigm(i).Outputs(3,1:5e4) = 0;
    
    ControlParadigm(i).Outputs(4,end-5e4:end) = 0;         
	ControlParadigm(i).Outputs(3,end-5e4:end) = 0;

end


% bundle the script that made this ControlParadigm with the controlparadigm
genScript = fileread(strcat(mfilename,'.m'));

save('MeanSwitching_2MFC_Kontroller_paradigm.mat','ControlParadigm','genScript');




