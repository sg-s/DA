% projectStimulusm.m
% method defined on the ORNdata class
% used to project the stimulus
% 
% created by Srinivas Gorur-Shandilya at 8:30 , 24 January 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.



function obj = projectStimulus(obj,projection_type)

assert(nargin==2,'Incorrect # of input arguments')
assert(ischar(projection_type),'2nd argument should be a string')

if strcmp(projection_type,'firing')
	obj.firing_projected = NaN*obj.stimulus;
	time = obj.dt*(1:length(obj.stimulus));
	for i = 1:obj.n_trials
	   obj.firing_projected(:,i) = convolve(time,obj.stimulus(:,i),obj.K_firing(:,i),obj.filtertime_firing);
	end
elseif strcmp(projection_type,'LFP')
	obj.LFP_projected = NaN*obj.stimulus;
	time = obj.dt*(1:length(obj.stimulus));
	for i = 1:obj.n_trials
	   obj.LFP_projected(:,i) = convolve(time,obj.stimulus(:,i),obj.K_LFP(:,i),obj.filtertime_LFP);
	end
else
	error('I dont know what you want me to project')
end 




