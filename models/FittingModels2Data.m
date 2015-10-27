% FittingModels2Data
% 
% created by Srinivas Gorur-Shandilya at 10:37 , 20 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

% add homebrew path
path1 = getenv('PATH');
if isempty(strfind(path1,':/usr/local/bin'))
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

load('/code/da/data/MSG_per_neuron')

% convert data into a nicer format
clear data
data = struct;
for i = 1:8
	data(i).stimulus = mean2([MSG_data(i,:).stim]);
	data(i).response = mean2([MSG_data(i,:).resp]);
end


%% Integral Feedback Model
% This model comes from the eLife paper on ORN responses in larvae, where the authors show that it performs well for very non-Gaussian odor stimuli (odor ramps, etc). This is how it looks:
% 
% <</code/da/models/elife-model.png>>
%


clear p
p.   b1= 4.3862;
p.   b2= 0.1238;
p.   b3= 9.9415e-04;
p.   b4= 0.2937;
p.   b5= 0.0087;
p.   a2= 0.0011;
p.theta= 0.6676;
p.   a3= 0.2984;

characteriseModel(@IFB,p,data)

prettyFig('fs=14;')

if being_published
	snapnow
	delete(gcf)
end


%% Incoherent Feed-Forward Model
% This also comes from the same paper.

clear p
p.   b1 = 9.9004;
p.   b2 = 6.5817;
p.   b3 = 0.0052;
p.   b4 = 0.2777;
p.   b5 = 0.0092;
p.   a1 = 5.4297;
p.   a2 = 0.0028;
p.theta = 0.4964;

characteriseModel(@IFF,p,data)

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
if being_published
	snapnow
	delete(gcf)
end



%% Methylation-Binding Model
% In this model, we consider receptors that can be bound to receptors or unbound, and can be methylated or unmethylated. We then write an abstract representation of the fraction of receptors bound and methylated:
%
% <</code/da/models/biophysical-model-v4.png>>
%

clear p
p.    sb =  1.0601e+04;
p.    sm =  585.1269;
p.    kb =  0.4404;
p.    km =  5.8281;
p.hill_A =  488.8477;
p.hill_k =  0.5152;
p. theta =  0.2584;

characteriseModel(@biophysicalModelv4,p,data)

prettyFig('fs=14;')

if being_published
	snapnow
	delete(gcf)
end

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

characteriseModel(@biophysicalModelv2,p,data)

prettyFig('fs=14;')

if being_published
	snapnow
	delete(gcf)
end





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
