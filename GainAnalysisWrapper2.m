% GainAnalysisWrapper.m
% A wrapper for Gain Analysis
% the purpose of this function is to do the gain analysis
% and make all the plots needed with just one line. no more fucking around.
% 
% created by Srinivas Gorur-Shandilya at 4:24 , 04 February 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function [p_LN,l,h,low_gof,high_gof,history_lengths] = GainAnalysisWrapper2(varargin)

if ~nargin
	help GainAnalysisWrapper
	return
elseif iseven(nargin)
    for ii = 1:nargin
        temp = varargin{ii};
        if ischar(temp)
            eval(strcat(temp,'=varargin{ii+1};'));
        end
    end
    clear ii
else
    error('Inputs need to be name value pairs')
end

% check for required inputs and prepare a structure to hash
x = struct;
if exist('response','var') && isvector(response)
	x.response = response;
else
	error('You need to specify a vector response')
end

if exist('stimulus','var') && isvector(stimulus)
	x.stimulus = stimulus;
else
	error('You need to specify a vector stimulus')
end

if exist('prediction','var') && isvector(prediction)
	x.prediction = prediction;
else
	error('You need to specify a vector prediction')
end

if exist('time','var') && isvector(time)
	x.time = time;
else
	error('You need to specify a vector time')
end

% check for optional arguments

if exist('history_lengths','var') 
else
	% find the correlation time of the stimulus
	[~,~,~,ct]=FindCorrelationTime(stimulus);
	dt = mean(diff(time));
	ct = ct*dt;

	% make history lengths based on the correlation time of the data
	history_lengths = (3*floor(1000*logspace(log10(ct/2),1,30)/3))/1e3;
	if nargin < 5 || isempty(example_history_length)
		example_history_length = history_lengths(10);
	end
end
x.history_lengths = history_lengths;

if exist('example_history_length','var') 
else
	example_history_length = history_lengths(10);
end

if exist('frac','var') 
else
	frac = .33;
end
x.frac = frac;

if exist('ph','var') 
else
	ph = [];
end

% clean up a little and ignore NaNs
x.filter_length = 499; % what does this even do??
rm_this = [find(isnan(response)) find(isnan(prediction)) ];
x.response(rm_this) = [];
x.prediction(rm_this) = [];
x.stimulus(rm_this) = [];
x.time(rm_this) = [];


% check cache to see if we have already computed this
hash = DataHash(x);
cached_data = cache(hash);
if isempty(cached_data)
	[p_LN,l,h,low_gof,high_gof] = GainAnalysis4(x,history_lengths,example_history_length,ph);
	cache(hash,p_LN);
	% also cache the example history length
	s=abs(l-h);
	s(p_LN(1,:)>0.05)=NaN;
	[~,loc]=max(s);
	ehl = history_lengths(loc);
	cache(DataHash(p_LN),ehl);

else
	% cached data exists. let's use that 
	p_LN = cached_data;
	if nargin < 5
		ehl = cache(DataHash(p_LN));
	else
		ehl = example_history_length;
	end
	[p_LN,l,h,low_gof,high_gof] = GainAnalysis4(x,history_lengths,ehl,ph,p_LN);
end

if ~isempty(ph)
	xlabel(ph(3),'Prediction (Hz)')
	ylabel(ph(3),'Data (Hz)')
	set(ph(4),'XScale','log')
end
