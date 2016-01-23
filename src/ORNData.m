% class definition for ORN data
% 
% created by Srinivas Gorur-Shandilya at 4:27 , 23 January 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

classdef ORNData
	properties
		firing_rate
		n_trials
		stimulus
		projected_stimulus
		K
		filtertime
		filter_length = 700;
		regularisation_factor = .1;
		dt = 1e-3;
		K_hash = 'filter not computed';
	end


	methods
   		function obj = set.firing_rate(obj,value)
   			% validate input
   			assert(isnumeric(value),'Firing rate should be numeric')

   			% if it is a matrix, it should be oriented properly
   			[value] = orientMatrix(value);

   			% it should match stimulus, if it exists
   			if ~isempty(obj.n_trials)
   				assert(obj.n_trials == size(value,2),'# of trials of firing rate should match # of trials of stimulus')
   			else
   				obj.n_trials = size(value,2);
   			end
   			obj.firing_rate = value;

   		end

   		function obj = set.stimulus(obj,value)
   			% validate input
   			assert(isnumeric(value),'Stimulus should be numeric')

   			% if it is a matrix, it should be oriented properly
   			[value] = orientMatrix(value);

   			% it should match stimulus, if it exists
   			if ~isempty(obj.n_trials)
   				assert(obj.n_trials == size(value,2),'# of trials of firing rate should match # of trials of stimulus')
   			else
   				obj.n_trials = size(value,2);
   			end
   			obj.stimulus = value;
   		end

   		function obj = set.K(obj,value)
   			% validate input
   			assert(isnumeric(value),'Filter should be numeric')

   			% if it is a matrix, it should be oriented properly
   			[value] = orientMatrix(value);
   			assert(size(value,2) == obj.n_trials,'Filter dimensions do not match # of trials')

   			obj.K = value;
   			%obj = projectStimulus(obj);
   		end

   		function obj = projectStimulus(obj)
   			obj.projected_stimulus = NaN*obj.stimulus;
   			time = obj.dt*(1:length(obj.stimulus));
   			for i = 1:obj.n_trials
   				obj.projected_stimulus(:,i) = convolve(time,obj.stimulus(:,i),obj.K(:,i),obj.filtertime);
   			end
   		end

   	   	function obj = extractFilter(obj)
   			if ~isempty(obj.firing_rate) && ~isempty(obj.stimulus) 
   				% use hashes to check if anything has changed since last compute
   				K_hash = dataHash([obj.firing_rate(:); obj.stimulus(:); obj.regularisation_factor(:); obj.filter_length(:)]);
   				if ~strcmp(K_hash,obj.K_hash)
   					disp('computing filter...')
	   				K = zeros(obj.filter_length,obj.n_trials);
	   				for i = 1:obj.n_trials
	   					textbar(i,obj.n_trials)
	   					[temp, filtertime] = fitFilter2Data(obj.stimulus(:,i), obj.firing_rate(:,i),'filter_length',obj.filter_length+200,'reg',obj.regularisation_factor,'offset',100+100);
	   					K(:,i) = temp(101:end-100);
	   					filtertime = filtertime(101:end-100);
	   				end
	   				obj.K_hash = K_hash;
	   				obj.K = K;
	   				obj.filtertime = filtertime*obj.dt;
	   				obj = projectStimulus(obj);
	   			end
   			end
   		end

   		function [plot_handles] = plot(varargin)
   			temp = varargin{1};
   			if strcmp(class(temp),'matlab.graphics.axis.Axes')
   				plot_here = varargin{1};
   				varargin{1} = [];
   			else
   				plot_here = gca;
   			end
   			assert(length(varargin)>1,'Not enough input arguments.')
   			obj = varargin{1};
   			switch varargin{2}
   			case 'K' 
   				plot_handles = plot(plot_here,obj.filtertime,obj.K);
   				xlabel('Filter Lag (s)')
   				ylabel('Filter')
   			case 'firing_rate'
   			end 
   			
   		end
   	end	
end

function [m] = orientMatrix(m)
	sz = size(m);
	if sz(2) > sz(1)
		m = m';
	end

end