function [cost]  = XJW_cost_function(x,data,IgnoreInitial,RemoveMean)
% convert the inputs into the parameter array that XJWNeuronWrapper needs
if ~isstruct(x)
	p= ValidateXJWParameters(x);
end

% unpack data
stimulus = data.stimulus;
response = data.response;
time = data.time;

if nargin < 3
	IgnoreInitial = floor(length(data.stimulus)/20);
end

if nargin < 4
	RemoveMean = 0;
end


if isvector(stimulus)
	% just a 1D fitting
	% now find the result from the guess
	[~,Rguess] = XJWNeuronWrapper(time,stimulus,p);
	Rguess = Rguess(:);
	response = response(:);
	if RemoveMean
		Rguess = Rguess- mean(Rguess(IgnoreInitial:end));
	end

	cost = Cost2(response(IgnoreInitial:end),Rguess(IgnoreInitial:end));
else
	% we're simultaneously fitting many different things
	error('WARNIGN!!! NEED TO CODE ONLY_THESE_POINTS TO HANDLE DISPARATE DATA SETS')
	Rguess = 0*response;
	for i = 1:size(stimulus,2)
		Rguess(:,i) = DA_integrate2(stimulus(:,i),p);
	end
	clear i
	Rguess = Rguess(:);
	response = response(:);
	cost = CostFunctionHandle(response,Rguess);
	
end
