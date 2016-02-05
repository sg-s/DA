% fitParametricFilters.m
% extracts parametric filters for class ORN data
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 05 February 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [obj] = fitParametricFilters(obj)

% respect constrains on where to calculate this 
if isempty(obj.use_these_trials)
	utt = 1:obj.n_trials;
else
	utt = obj.use_these_trials;		
end
if isempty(obj.use_this_segment)
	uts = 1:length(obj.stimulus);
else
	uts = obj.use_this_segment;
end	


% do the firing rate
if ~isempty(obj.firing_rate) && ~isempty(obj.stimulus) 

	assert(min(obj.filtertime_firing) >= 0, 'You need to specify a all positive filtertime vector')

	% use hashes to check if anything has changed since last compute
	K_firing_hash = dataHash([vectorise(obj.firing_rate(uts,utt)); vectorise(obj.stimulus(uts,utt)); obj.filtertime_firing(:)]);
	if ~strcmp(K_firing_hash,obj.K_firing_hash)
		filter_length = length(obj.filtertime_firing);
		clear p
		p.     n = 0.2188;
		p.     A = 0.5977;
		p.  tau1 = 186.7812;
		p.  tau2 = 114.6250;
		p. scale = 0.4618;
		p.offset = 11.0639;
		K = NaN(filter_length,obj.n_trials);
		for i = 1:obj.n_trials
			if  ismember(i,utt)
				clear d
				d.stimulus = obj.stimulus(uts,i); 
				d.stimulus = d.stimulus - mean(d.stimulus);  % Baccus-Meister norm.
				d.stimulus = d.stimulus/std(d.stimulus);
				d.response = obj.firing_rate(uts,i);
				d.response = d.response - mean(d.response);  % Baccus-Meister norm.
				d.response = d.response/std(d.response);
				p = fitModel2Data(@pLinearModel,d,'nsteps',500,'p0',p);
				K(:,i) = filter_gamma2(obj.filtertime_firing/obj.dt,p);
				K(:,i) = K(:,i)*p.scale;

			end
		end
		obj.K_firing_hash = K_firing_hash;
		obj.K_firing = K;
		obj.regularisation_factor = NaN;
		obj = projectStimulus(obj,'firing');
	end
end