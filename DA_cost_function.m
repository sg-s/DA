% created by Srinivas Gorur-Shandilya at 12:27 , 05 December 2013. Contact me at http://srinivas.gs/contact/
% DA_cost_function is a function that evaluates the response predicted by the DA model
% to the stimulus and and compares it to the actual response. it calcualtes the absolute
% error of the prediction, so has a minimum of 0 when the prediction is perfect. 
function [cost]  = DA_cost_function(x,data,CostFunctionHandle,algo)
% convert the inputs into the parameter array that DA_integrate needs
if ~isstruct(x)
	p= ValidateDAParameters(x,algo);
else
	p = x;
end


stimulus = data.stimulus;
response = data.response;



if isvector(stimulus)
	% just a 1D fitting
	% now find the result from the guess
	Rguess = DA_integrate(stimulus,p);
	Rguess = Rguess(:);
	response = response(:);
	stimulus = stimulus(:);
	Rguess = Rguess - mean(Rguess(500:end));
	error('need some careful thought about how you calcualte the mean: throw out some initial segment?')

	cost = CostFunctionHandle(response,Rguess);
else
	% we're simultaneously fitting many different things
	Rguess = 0*response;
	for i = 1:size(stimulus,1)
		Rguess(i,:) = DA_integrate(stimulus(i,:),p);
	end
	clear i
	Rguess = Rguess';
	response = response';
	Rguess = Rguess(:);
	response = response(:);
	cost = CostFunctionHandle(response,Rguess);
	
end



