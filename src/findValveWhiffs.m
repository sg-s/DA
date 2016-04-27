%% findValveWhiffs.m
% finds intervals in a timeseries in a ORNData object when the valve just opened
% needs information about valve signal to work
% assumes that a valve opens and closes to deliver odour
% usage:
% [ons,offs] = findValveWhiffs(orn_data)
%
function [valve_ons,valve_offs] = findValveWhiffs(varargin)


% options and defaults
options.t_after_valve_opens = 100; % milliseconds

if nargout && ~nargin
	varargout{1} = options;
end

% grab the ORNData object from the arguments
for i = 1:length(varargin)
	if isa((varargin{i}),'ORNData')
		o = varargin{i};
		varargin(i) = [];
		break
	end
end

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

% find all valve ons and offs
[valve_ons,valve_offs] = computeOnsOffs(o.valve(o.use_this_segment));

% remove the last one because it might be at the end
valve_ons(end) = []; valve_offs(end) = [];

% make sure whiffs do not exceed a maximum length
valve_offs((valve_offs - valve_ons) > options.t_after_valve_opens) = valve_ons((valve_offs - valve_ons) > options.t_after_valve_opens) + options.t_after_valve_opens;

valve_ons = valve_ons + find(o.use_this_segment,1,'first');
valve_offs = valve_offs + find(o.use_this_segment,1,'first');