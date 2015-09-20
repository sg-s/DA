% GainAnalysisWrapper.m
% A wrapper for Gain Analysis
% the purpose of this function is to do the gain analysis
% and make all the plots needed with just one line. no more fucking around.
% 
% created by Srinivas Gorur-Shandilya at 4:24 , 04 February 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function [p_LN,l,h,low_gof,high_gof,history_lengths,handles] = GainAnalysisWrapper2(varargin)

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
if exist('verbosity','var') 
else
	verbosity = 0;
end

if exist('history_lengths','var') 
else
	% find the correlation time of the stimulus
	[~,~,~,ct]=FindCorrelationTime(stimulus);
	dt = mean(diff(time));
	ct = ct*dt;

	% make history lengths based on the correlation time of the data
	history_lengths = (3*floor(1000*logspace(log10(ct/2),1,30)/3))/1e3;
	if isempty(example_history_length)
		example_history_length = history_lengths(10);
	end
end
x.history_lengths = history_lengths;


if exist('frac','var') 
else
	frac = .33;
end
x.frac = frac;

if exist('ph','var') 
else
	ph = [];
end

if exist('use_cache','var') 
else
	use_cache = true;
end

% new -- now allows you to choose which engine to use. 
if exist('engine','var') 
else
	engine = @GainAnalysis4;
end


% clean up a little and ignore NaNs
rm_this = [find(isnan(response)); find(isnan(prediction)) ];

assert(length(rm_this)<length(response)/2,'GainAnalysisWrapper2: It looks like most of the response/prediction are NaN.')

x.response(rm_this) = [];
x.prediction(rm_this) = [];
x.stimulus(rm_this) = [];
x.time(rm_this) = [];
x.engine = engine;

% check cache to see if we have already computed this
hash = dataHash(x);
cached_data = cache(hash);
if ~use_cache
	cached_data = [];
end
if isempty(cached_data)
	if exist('example_history_length','var') 
	else
		if verbosity
			disp('using default example histroy length')
		end
		example_history_length = history_lengths(10);
	end
	[p_LN,l,h,low_gof,high_gof,~,~,handles] = x.engine(x,history_lengths,example_history_length,ph);
	cache(hash,p_LN);
	% also cache the example history length
	g = l-h;
	g(~p_LN(2,:)) = -Inf;
	[~,loc]=max(g);
	ehl = history_lengths(loc);

	cache(dataHash(p_LN),ehl);
	if verbosity
		disp('Computing the best example history length for next time. and that is...')
		disp(ehl)
	end

else
	if verbosity
		disp('cached data exists. lets use that ')
	end
	p_LN = cached_data;
	if exist('example_history_length','var') 
		ehl = example_history_length;
	else
		ehl = cache(dataHash(p_LN));
	end
	[p_LN,l,h,low_gof,high_gof,~,~,handles] = x.engine(x,history_lengths,ehl,ph,p_LN);
end

try
	xlabel(ph(3),'Prediction (Hz)')
	ylabel(ph(3),'Data (Hz)')
catch
end
try
	set(ph(4),'XScale','log')
catch
end

