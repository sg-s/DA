% STA.m
% computes the spike triggered averaged stimulus.
% works with one or many trials
% usage:
% K = STA(spikes,stimulus,before,after)
% 
% where 
% 
% spikes is a sparse matrix with 1 indicating a spike
% stimulus is a matrix the same size as spikes
% before is a positive integer 
% after is a positive integer
% K is a matrix that is before+1+after elements long
% 
% created by Srinivas Gorur-Shandilya at 3:48 , 18 March 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function K = STA(spikes,stimulus,before,after,normalise)

if ~nargin
	help STA
	return
elseif nargin < 2
	error('Not enough input arguments')
elseif nargin < 3
	before = 10000;
	after = 6000; 
elseif nargin < 5
	normalise= true;
end

if ~issparse(spikes)
	error('First argument should be a sparse matrix')
end

if length(spikes) ~= length(stimulus)
	error('first two arguments have to be the same length')
end

% orient matrices properly
if ~isvector(spikes)
	if size(spikes,2) > size(spikes,1)
		spikes = spikes';
	end
	
else
	spikes = spikes(:);
end
if ~isvector(stimulus)
	if size(stimulus,2) > size(stimulus,1)
		stimulus = stimulus';
	end
else
	stimulus = stimulus(:);
end
K = zeros(before+after+1,width(spikes));

for i = 1:width(spikes)
	% normalise if needed
	if normalise
		stimulus(:,i) =  stimulus(:,i) - mean(stimulus(:,i));
		stimulus(:,i) =  stimulus(:,i)/std(stimulus(:,i));
	end
	permitted_spikes = find(spikes(:,i));
	permitted_spikes(permitted_spikes<before) = [];
	permitted_spikes(permitted_spikes>length(stimulus)-after) = [];
	for j = 1:length(permitted_spikes)
		this_spike = permitted_spikes(j);
		K(:,i) = K(:,i) + stimulus(this_spike-before:this_spike+after,i);
	

		
	end
	K(:,i) = K(:,i)/length(permitted_spikes);

end



