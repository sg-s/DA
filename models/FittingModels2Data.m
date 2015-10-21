% FittingModels2Data
% 
% created by Srinivas Gorur-Shandilya at 10:37 , 20 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

% add homebrew path
path1 = getenv('PATH');
if isempty(strfind(path1,[':/usr/local/bin']))
    path1 = [path1 ':/usr/local/bin'];
end
setenv('PATH', path1);

% this code determines if this function is being called
calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
		unix(['tag -a publish-failed ',which(mfilename)]);
		unix(['tag -r published ',which(mfilename)]);
	end
end
tic

%% Fitting Models to Data
% In this document we fit some models to data we have, and study if these models can account for gain and kinetics we observe in the data. 

%% Simple Biophysical Model 
% The first model we will do is a simple biophysical model involving receptor binding, a slow diffusible factor, and a static output nonlinearity:
% 
% <</code/da/models/biophysical_model.png>>
%


load('/code/da/data/MSG_per_neuron')

clear p
p.    kb = 102.8451;
p.    sb = 22.8818;
p.    ko = 126.4313;
p.    kc = 231.5404;
p.     a = 228.7689;
p.     b = 1000;
p.hill_A = 44.7539;
p.hill_n = 6.2676;
p.hill_k = 0.3300;
p.     n = 2.3000;

% convert data into a nicer format
clear data
for i = 1:8
	data(i).stimulus = mean2([MSG_data(i,:).stim]);
	data(i).response = mean2([MSG_data(i,:).resp]);
end

characteriseModel(@biophysicalModelv2,p,data)

prettyFig('fs=14;')

if being_published
	snapnow
	delete(gcf)
end

%% Dynamical Adaptation Model
% In this section we do the same analysis for a reduced version of the DA model, which consists of two filters, one dividing the other. 


clear p
p.    s0= -0.3502;
p.   n_z= 2;
p. tau_z= 165.5000;
p.   n_y= 2;
p. tau_y= 27.0156;
p.     C= 0.2001;
p.     A= 357.2305;
p.     B= 14.4248;

characteriseModel(@DAModelv2,p,data)

prettyFig('fs=14;')


%% Version Info
% The file that generated this document is called:
disp(mfilename)

%%
% and its md5 hash is:
Opt.Input = 'file';
disp(dataHash(strcat(mfilename,'.m'),Opt))

%%
% This file should be in this commit:
[status,m]=unix('git rev-parse HEAD');
if ~status
	disp(m)
end

t = toc;

%% 
% This document was built in: 
disp(strcat(oval(t,3),' seconds.'))

% tag the file as being published 

if being_published
	unix(['tag -a published ',which(mfilename)]);
	unix(['tag -r publish-failed ',which(mfilename)]);
end
