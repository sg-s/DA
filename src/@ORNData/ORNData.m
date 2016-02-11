% class definition for ORN data
% 
% created by Srinivas Gorur-Shandilya at 4:27 , 23 January 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

classdef ORNData
	properties
      % raw source data
		firing_rate           % firing rate (in Hz. a matrix)
		stimulus              % same dimensions as firing_rate (in V)
      LFP
      n_trials              % # of trials. integer. 
      valve
      MFC_control
      MFC_signal
      spikes

      % timing information
      dt = 1e-3;
      timescale_inst_gain = 50; % in ms
      filtertime_firing = 1e-3*(1:700) - 100*1e-3;
      filtertime_LFP = 1e-3*(1:700) - 100*1e-3;

      % filters and projections
      regularisation_factor = .1;
      K_firing
      K_LFP
      K_firing_hash = 'filter not computed';
      K_LFP_hash = 'filter not computed';
      firing_projected
      LFP_projected

      % gain
      gain_LFP 
      gain_firing

      inst_gain_LFP 
      inst_gain_LFP_err          
      inst_gain_firing     % averaged over all trials
      inst_gain_firing_err
      

      % metadata
      original_name
      data_creator
      neuron_name
      odour_name
      orn % ID of orn
      fly % ID of fly
      paradigm 
      console_log = '';

      % filters to plot subparts of the data
      use_these_trials
      use_this_segment






	end


	methods
      % only set/get methods are here; other methods in their own m files

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
         obj.console_log = [obj.console_log char(10) ' ' datestr(now) '    Firing Rate set'];


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
         obj.console_log = [obj.console_log char(10) ' ' datestr(now) '   Stimulus set'];
         

   	end

   	function obj = set.K_firing(obj,value)
   		% validate input
   		assert(isnumeric(value),'Filter should be numeric')

   		% if it is a matrix, it should be oriented properly
   		[value] = orientMatrix(value);
   		assert(size(value,2) == obj.n_trials,'Filter dimensions do not match # of trials')

   		obj.K_firing = value;
         obj.console_log = [obj.console_log char(10) ' ' datestr(now) '   K_firing rate set'];

   		obj = projectStimulus(obj,'firing');
   	end

      function obj = set.regularisation_factor(obj,value)
         % validate input
         assert(isnumeric(value),'regularisation_factor must be a +ve number')
         assert(isscalar(value),'regularisation_factor must be a +ve number')
         assert(value>0,'regularisation_factor must be a +ve number')
         obj.regularisation_factor = value;
         %obj = backOutFilters(obj);
      end


      function obj = set.filtertime_LFP(obj,value)
         % validate input
         assert(isvector(value),'filtertime_LFP must be a vector')

         obj.filtertime_LFP = value;
         %obj = backOutFilters(obj);
      end

      function obj = set.use_this_segment(obj,value)
         % validate input
         assert(length(value) == length(obj.stimulus),'use_this_segment must be a vector that is as long as the stimulus')
         obj.use_this_segment = logical(value);
      end

   end

end


function [m] = orientMatrix(m)
   sz = size(m);
   if sz(2) > sz(1)
      m = m';
   end
end