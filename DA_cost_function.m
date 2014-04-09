% created by Srinivas Gorur-Shandilya at 12:27 , 05 December 2013. Contact me at http://srinivas.gs/contact/
% DA_cost_function is a function that evaluates the response predicted by the DA model
% to the stimulus and and compares it to the actual response. it calcualtes the absolute
% error of the prediction, so has a minimum of 0 when the prediction is perfect. 
function [cost]  = DA_cost_function(x,PID,f,CostFunctionHandle)
% convert the inputs into the parameter array that DA_integrate needs
p= ValidateDAParameters(x);

% now find the result from the guess
Rguess = DA_integrate(PID,p);

f = f(:);
Rguess = Rguess(:);


cost = CostFunctionHandle(f,Rguess);
