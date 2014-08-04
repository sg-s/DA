% Fit NonLinear Function to Model
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function [NFit,x]  = FitNModel(data,x0)

switch nargin
case 0
	help FitNModel
	return
case 1
	lb = [0 0 0];
	ub = [max(data.response) 10 10];
	x0 = [max(data.response) 1 1];
case 2
	lb = x0/10;
	ub = x0*10;
end

% The Hill function is not defined for negative values, so just remove the lowest
data.stimulus = data.stimulus - min(data.stimulus);

psoptions = psoptimset('UseParallel',true,'CompletePoll', 'on', 'Vectorized', 'off','Display','iter','MaxIter',100,'MaxFunEvals',10000);
x = patternsearch(@(x) NCostFunction(x,data),x0,[],[],[],[],lb,ub,psoptions);	


% get the final output
NFit = hill(x,data.stimulus);

% debug
figure, hold on
plot(data.response), hold on
plot(NFit,'r')

disp(rsquare(data.response,NFit))