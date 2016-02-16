% ORNData/fitParametricFilters.m
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

	% use hashes to check if anything has changed since last compute
	K_firing_hash = dataHash([vectorise(obj.firing_rate(uts,utt)); vectorise(obj.stimulus(uts,utt)); obj.filtertime_firing(:)]);
	if ~strcmp(K_firing_hash,obj.K_firing_hash)
		filter_length = length(obj.filtertime_firing);
		K = NaN(filter_length,obj.n_trials);
		for i = 1:obj.n_trials
			if  ismember(i,utt)
				clear d p
				p.     n = 1;
				p.     A = .5;
				p.  tau1 = 20;
				p.  tau2 = 100;
				p. scale = 0.04;
				p.offset = 0;
				p.x_offset = -10;
				d.stimulus = obj.stimulus(uts,i); 
				d.stimulus = d.stimulus - mean(d.stimulus);  % Baccus-Meister norm.
				d.stimulus = d.stimulus/std(d.stimulus);
				d.response = obj.firing_rate(uts,i);
				d.response = d.response - mean(d.response);  % Baccus-Meister norm.
				d.response = d.response/std(d.response);
				p = fitModel2Data(@pLinearModel,d,'nsteps',100,'p0',p);		
				[~,temp] = pLinearModel(d.stimulus,p);
				if p.x_offset == 0
					disp('no offset')
					keyboard
					K(:,i) = interp1(1:length(temp),temp,obj.filtertime_firing,'nearest','extrap');

				elseif p.x_offset > 0
					temp = [zeros(round(p.x_offset),1); temp(:)];
					ft = ((1:length(temp)) - p.x_offset)*obj.dt;
					K(:,i) = interp1(ft,temp,obj.filtertime_firing,'nearest','extrap');
				else
					ft = ((1:length(temp)) + p.x_offset)*obj.dt;
					K(:,i) = interp1(ft,temp,obj.filtertime_firing,'nearest','extrap');
				end

				% get the scale right
				K(:,i) = K(:,i)*p.scale;
				

			end
		end
		obj.K_firing_hash = K_firing_hash;
		obj.K_firing = K;
		obj.regularisation_factor = Inf;
		obj = projectStimulus(obj,'firing');
	end
end

% do the LFP
if ~isempty(obj.LFP) && ~isempty(obj.stimulus) 

	% use hashes to check if anything has changed since last compute
	K_LFP_hash = dataHash([vectorise(obj.LFP(uts,utt)); vectorise(obj.stimulus(uts,utt)); obj.filtertime_LFP(:)]);
	if ~strcmp(K_LFP_hash,obj.K_LFP_hash)
		filter_length = length(obj.filtertime_LFP);
		K = NaN(filter_length,obj.n_trials);
		for i = 1:obj.n_trials
			if  ismember(i,utt)
				
				% fit pLFP.m
				clear p d
				p.tau = 10; p.A = .1; p.offset = 0; p.x_offset = 0; p.n = 2; 
				d.stimulus = obj.stimulus(uts,i); 
				d.stimulus = d.stimulus - mean(d.stimulus);  % Baccus-Meister norm.
				d.stimulus = d.stimulus/std(d.stimulus);
				d.response = obj.LFP(uts,i);
				d.response = d.response - mean(d.response);  % Baccus-Meister norm.
				d.response = d.response/std(d.response);
				p = fitModel2Data(@pLFP,d,'nsteps',100,'p0',p);		
				[~,temp] = pLFP(d.stimulus,p);
				if p.x_offset == 0
					disp('zero offset for LFP')
					keyboard
					K(:,i) = interp1(1:length(temp),temp,obj.filtertime_LFP,'nearest','extrap');
				elseif p.x_offset > 0
					temp = [zeros(round(p.x_offset),1); temp(:)];
					ft = ((1:length(temp)) - p.x_offset)*obj.dt;
					K(:,i) = interp1(ft,temp,obj.filtertime_LFP,'nearest','extrap');
				else
					ft = ((1:length(temp)) + p.x_offset)*obj.dt;
					K(:,i) = interp1(ft,temp,obj.filtertime_LFP,'nearest','extrap');
				end

				% get the scale right
				K(:,i) = K(:,i)*-p.A;
				

			end
		end
		obj.K_LFP_hash = K_LFP_hash;
		obj.K_LFP = K;
		obj.regularisation_factor = Inf;
		obj = projectStimulus(obj,'firing');
	end
end