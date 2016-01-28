% extracts filters for class ORNData
% 
% created by Srinivas Gorur-Shandilya at 6:30 , 23 January 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function obj = backOutFilters(obj,filter_type)

if nargin < 2
	filter_type = 'all';
end

assert(nargin<3,'Too many input arguments')
assert(ischar(filter_type),'2nd argument should be a string: {"all","LFP","firing"}');

if  strcmp(filter_type,'all') || strcmp(filter_type,'firing')
	% do the firing rate
	if ~isempty(obj.firing_rate) && ~isempty(obj.stimulus) 
		% use hashes to check if anything has changed since last compute
		K_firing_hash = dataHash([obj.firing_rate(:); obj.stimulus(:); obj.regularisation_factor(:); obj.filtertime_firing(:)]);
		if ~strcmp(K_firing_hash,obj.K_firing_hash)
			disp('computing filter...')
			K = zeros(obj.filter_length,obj.n_trials);
			for i = 1:obj.n_trials
				textbar(i,obj.n_trials)
				[temp, filtertime] = fitFilter2Data(obj.stimulus(:,i), obj.firing_rate(:,i),'filter_length',obj.filter_length+200,'reg',obj.regularisation_factor,'offset',obj.filter_offset+100);
				K(:,i) = temp(101:end-100);
				filtertime = filtertime(101:end-100);
			end
			obj.K_firing_hash = K_firing_hash;
			obj.K_firing = K;
			obj.filtertime_firing = filtertime*obj.dt;
			obj = projectStimulus(obj,'firing');
		end
	end
end
if strcmp(filter_type,'all') || strcmp(filter_type,'LFP')
	% do the LFP
	if ~isempty(obj.LFP) && ~isempty(obj.stimulus) 
		% use hashes to check if anything has changed since last compute
		K_LFP_hash = dataHash([obj.LFP(:); obj.stimulus(:); obj.regularisation_factor(:); obj.filtertime_LFP(:)]);
		if ~strcmp(K_LFP_hash,obj.K_LFP_hash)
			disp('computing filter...')
			filter_length = length(obj.filtertime_LFP);
			filter_offset = find(obj.filtertime_LFP == 0);
			K = zeros(filter_length,obj.n_trials);
			for i = 1:obj.n_trials
				textbar(i,obj.n_trials)
				[temp, filtertime] = fitFilter2Data(obj.stimulus(:,i), obj.LFP(:,i),'filter_length',filter_length+200,'reg',obj.regularisation_factor,'offset',filter_offset+100);
				K(:,i) = temp(101:end-100);
				filtertime = filtertime(101:end-100);
			end
			obj.K_LFP_hash = K_LFP_hash;
			obj.K_LFP = K;
			obj.filtertime_firing = filtertime*obj.dt;
			obj = projectStimulus(obj,'LFP');
		end
	end
else
	error('I cant understand what you want me to back out')
end