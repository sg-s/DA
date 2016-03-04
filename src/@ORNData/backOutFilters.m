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

% recursively call this to work on arrays of objects
original_dimensions = size(obj);
obj = obj(:);
if length(obj) > 1
	for i = 1:length(obj)
		obj(i) = backOutFilters(obj(i),filter_type);
	end
	obj = reshape(obj,original_dimensions);
	return
end

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


assert(nargin<3,'Too many input arguments')
assert(ischar(filter_type),'2nd argument should be a string: {"all","LFP","firing"}');

if  strcmp(filter_type,'all') || strcmp(filter_type,'firing')
	% do the firing rate
	if ~isempty(obj.firing_rate) && ~isempty(obj.stimulus) 
		% use hashes to check if anything has changed since last compute
		K_firing_hash = dataHash([vectorise(obj.firing_rate(uts,utt)); vectorise(obj.stimulus(uts,utt)); obj.regularisation_factor(:); obj.filtertime_firing(:)]);
		if ~strcmp(K_firing_hash,obj.K_firing_hash)
			filter_length = length(obj.filtertime_firing);
			filter_offset = find(obj.filtertime_firing == 0);
			disp('computing filter...')
			K = NaN(filter_length,obj.n_trials);

			if any(utt)
				for i = 1:obj.n_trials
					if  ismember(i,find(utt))
						[temp, filtertime] = fitFilter2Data(obj.stimulus(uts,i), obj.firing_rate(uts,i),'filter_length',filter_length+200,'reg',obj.regularisation_factor,'offset',filter_offset+100);
						K(:,i) = temp(101:end-100);
						filtertime = filtertime(101:end-100);
					end
				end
				obj.K_firing_hash = K_firing_hash;
				obj.K_firing = K;
				obj.filtertime_firing = filtertime*obj.dt;
				obj = projectStimulus(obj,'firing');
			end
		end
	end
end
if strcmp(filter_type,'all') || strcmp(filter_type,'LFP')
	% do the LFP
	if ~isempty(obj.LFP) && ~isempty(obj.stimulus) 
		% use hashes to check if anything has changed since last compute
		K_LFP_hash = dataHash([vectorise(obj.LFP(uts,utt)); vectorise(obj.stimulus(uts,utt)); obj.regularisation_factor(:); obj.filtertime_firing(:)]);
		if ~strcmp(K_LFP_hash,obj.K_LFP_hash)
			disp('computing filter...')
			filter_length = length(obj.filtertime_LFP);
			filter_offset = find(obj.filtertime_LFP == 0);
			K = NaN(filter_length,obj.n_trials);
			if any(utt)
				for i = 1:obj.n_trials
					if  ismember(i,find(utt))
						[temp, filtertime] = fitFilter2Data(obj.stimulus(uts,i), obj.LFP(uts,i),'filter_length',filter_length+200,'reg',obj.regularisation_factor,'offset',filter_offset+100);
						K(:,i) = temp(101:end-100);
						filtertime = filtertime(101:end-100);
					end
				end
				obj.K_LFP_hash = K_LFP_hash;
				obj.K_LFP = K;
				obj.filtertime_firing = filtertime*obj.dt;
				obj = projectStimulus(obj,'LFP');
			end
		end
	end
else
	error('I cant understand what you want me to back out')
end