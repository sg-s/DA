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
options.min_inst_gain_firing = 10;

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

% grab the response
resp = nanmean(o.firing_rate,2);
resp = resp(o.use_this_segment);

% find all valve ons and offs
[valve_ons,valve_offs] = computeOnsOffs(o.valve(o.use_this_segment));

% remove the last one because it might be at the end
valve_ons(end) = []; valve_offs(end) = [];
resp(resp<options.min_inst_gain_firing) = NaN;



for i = 1:length(valve_ons)
	a = find(~isnan(resp(valve_ons(i):valve_offs(i))),1,'first');
	if ~isempty(a)
		a = a + valve_ons(i);
		[~,z] = max(resp(a:valve_offs(i))); z = z + a;
		if z - a > options.t_after_valve_opens
			z = a + options.t_after_valve_opens;
		end
		if any(isnan(resp(a:z)))
			z = a+find(~isnan(resp(a:z)),1,'last');
		end
		valve_ons(i) = a;
		try
			valve_offs(i) = z;
		catch
			valve_offs(i) = NaN;
			valve_ons(i) = NaN;
		end

	end
end
valve_offs(isnan(valve_offs)) = [];
valve_ons(isnan(valve_ons)) = [];