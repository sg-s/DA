% created by Srinivas Gorur-Shandilya at 12:27 , 05 December 2013. Contact me at http://srinivas.gs/contact/
% DA_cost_function is a function that evaluates the response predicted by the DA model
% to the stimulus and and compares it to the actual response. it calcualtes the absolute
% error of the prediction, so has a minimum of 0 when the prediction is perfect. 
function [cost]  = DA_cost_function(x,data,CostFunctionHandle)
% convert the inputs into the parameter array that DA_integrate needs
if ~isstruct(x)
	p= ValidateDAParameters2(x);
end

% unpack data
stimulus = data.stimulus;
response = data.response;



if isvector(stimulus)
	% just a 1D fitting
	% now find the result from the guess
	Rguess = DA_integrate2(stimulus,p);
	Rguess = Rguess(:);
	response = response(:);
	stimulus = stimulus(:);
	
	cost = CostFunctionHandle(response,Rguess);
else
	% we're simultaneously fitting many different things
	% WARNIGN!!! NEED TO CODE ONLY_THESE_POINTS TO HANDLE DISPARATE DATA SETS
	Rguess = 0*response;
	for i = 1:size(stimulus,2)
		Rguess(:,i) = DA_integrate2(stimulus(:,i),p);
	end
	clear i
	Rguess = Rguess(:);
	response = response(:);
	cost = CostFunctionHandle(response,Rguess);
	
end



