% STA.m
% computes the spike triggered averaged stimulus.
% works with one or many trials
% minimal usage:
% K = STA(spikes,stimulus)
% 
% full usage:
% K = STA(spikes,stimulus,'normalise',true,'regulariseParameter',1,'before',100,'after',10);
% 
% created by Srinivas Gorur-Shandilya at 3:48 , 18 March 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function K = STA(spikes,stimulus,varargin)

% defaults
before = 1e4;
after = 1e3;
regulariseParameter = 0;
normalise = true;

if ~nargin
    help STA
    return
else
    if iseven(nargin)
    	for ii = 3:2:length(varargin)-1
        	temp = varargin{ii};
        	if ischar(temp)
            	eval(strcat(temp,'=varargin{ii+1};'));
        	end
    	end
	else
    	error('Inputs need to be name value pairs')
	end
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
K = K(1:10:end,:);

for i = 1:width(spikes)
	% normalise if needed
	if normalise
		stimulus(:,i) =  stimulus(:,i) - mean(stimulus(:,i));
		stimulus(:,i) =  stimulus(:,i)/std(stimulus(:,i));
	end
	permitted_spikes = find(spikes(:,i));
	permitted_spikes(permitted_spikes<before) = [];
	permitted_spikes(permitted_spikes>length(stimulus)-after) = [];

	Y = zeros(length(permitted_spikes),before+after+1);

	for j = 1:length(permitted_spikes)
		this_spike = permitted_spikes(j);
		Y(j,:) = stimulus(this_spike-before:this_spike+after,i);		
	end

	% downsample X
	Y = Y(:,1:10:end);

	if regulariseParameter ~= 0
		all_times = (floor(before)/10+1:10:length(stimulus)-floor(after/10));
		X = zeros(length(all_times),floor(before/10)+floor(after/10)+1);
		for j = 1:length(all_times)
			textbar(j,length(all_times))
			this_spike = all_times(j);
			X(j,:) = stimulus(this_spike-floor(before/10):this_spike+floor(after/10),i);		
		end

		C = X'*X;
		MeanEigenValue = trace(C)/length(C);

		K(:,i) = (mean(((C + regulariseParameter*MeanEigenValue*eye(length(C)))\Y')')); 
		K(:,i) = K(:,i)/max(K(:,i));
		K(:,i) = K(:,i)*max(mean2(Y));
	else
		K(:,i) = mean2(Y);
	end



end



