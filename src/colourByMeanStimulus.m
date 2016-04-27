%% colourByMeanStimulus.m
% builds a colour index based on the mean stimulus in some history window
% 
function [varargout] = colourByMeanStimulus(S,use_this_segment,varargin)

% options and defaults
options.history_length = 300;
options.n_colours = 100;

if nargout && ~nargin 
	varargout{1} = options;
    return
end

% validate and accept options
if iseven(length(varargin))
	for ii = 1:2:length(varargin)-1
	temp = varargin{ii};
    if ischar(temp)
    	if ~any(find(strcmp(temp,fieldnames(options))))
    		disp(['Unknown option: ' temp])
    		disp('The allowed options are:')
    		disp(fieldnames(options))
    		error('UNKNOWN OPTION')
    	else
    		options = setfield(options,temp,varargin{ii+1});
    	end
    end
end
elseif isstruct(varargin{1})
	% should be OK...
	options = varargin{1};
else
	error('Inputs need to be name value pairs')
end

% show the response vs. the projected stimulus
shat = computeSmoothedStimulus(S,options.history_length);
shat = shat-min(shat(use_this_segment));
shat = shat/max(shat(use_this_segment));
shat = 1+ceil(shat*options.n_colours-1);
shat(isnan(shat)) = 1;
shat(shat<1) = 1;
shat(shat>options.n_colours-1) = options.n_colours;

cc = parula(options.n_colours);
varargout{1} = cc(shat,:);