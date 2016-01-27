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
      filter_length = 700;
      filter_offset = 100;
      filtertime_firing = 1e-3*(1:700) - 100*1e-3;
      a = 10e3; % where do we start looking at the data
      z

      % filters, projections and gain
      regularisation_factor = .1;
      K_firing
      K_LFP
      K_firing_hash = 'filter not computed';
      K_LPF_hash = 'filter not computed';
      firing_projected
      LFP_projected
      inst_gain_LFP
      gain_LFP             % trial-wise 
      inst_gain_firing     % averaged over all trials
      inst_gain_firing_err
      gain_firing

      % metadata
      original_name
      data_creator
      neuron_name
      odour_name
      orn % ID of orn
      fly % ID of fly
      paradigm 






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

         % set z if unset 
         if isempty(obj.z)
            obj.z = length(value);
         end

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
         
         % set z if unset 
         if isempty(obj.z)
            obj.z = length(value);
         end
   	end

   	function obj = set.K_firing(obj,value)
   		% validate input
   		assert(isnumeric(value),'Filter should be numeric')

   		% if it is a matrix, it should be oriented properly
   		[value] = orientMatrix(value);
   		assert(size(value,2) == obj.n_trials,'Filter dimensions do not match # of trials')

   		obj.K_firing = value;
   		obj = projectStimulus(obj,'firing');
   	end

      function obj = set.regularisation_factor(obj,value)
         % validate input
         assert(isnumeric(value),'regularisation_factor must be a +ve number')
         assert(isscalar(value),'regularisation_factor must be a +ve number')
         assert(value>0,'regularisation_factor must be a +ve number')

         disp('Recomputing filters because you changed the regularisation_factor...')
         obj.regularisation_factor = value;
         obj = backOutFilters(obj);

      end

   end

end


function [m] = orientMatrix(m)
   sz = size(m);
   if sz(2) > sz(1)
      m = m';
   end
end