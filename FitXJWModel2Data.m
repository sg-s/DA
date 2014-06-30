% FitXJWModel2Data
% fits the Liu and Wang model (XJWNeuronWrapper) to ORN data. uses patternsearch and the optimisation toolbox. 
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [p,Rguess] = FitXJWModel2Data(data,x0,lb,ub,IgnoreInitial)


switch nargin 
	case 0
		help FitXJWModelToData
		return
	case 1
		x0 = [4.5  0.5   3.1   19    0.35   -58  -55 -57 -80  11899];
		lb = [1e-3 1e-4  0     1e-3  0      -80  -55 -99 -80   1 ];
		ub = [10   10    1e2   40    20     -40  -55 -55 -80  1e6];
		IgnoreInitial = 1;
	case 2
		lb = x0/2;
		ub = x0*10;
		IgnoreInitial = floor(length(data.stimulus)/20);
	case 3
		ub = x0*10;
		IgnoreInitial = floor(length(data.stimulus)/20);
	case 4
		IgnoreInitial = floor(length(data.stimulus)/20);
end



psoptions = psoptimset('UseParallel',true,'CompletePoll', 'on', 'Vectorized', 'off','Display','iter','MaxIter',100,'MaxFunEvals',10000);
x = patternsearch(@(x) XJW_cost_function(x,data,IgnoreInitial),x0,[],[],[],[],lb,ub,psoptions);
p = ValidateXJWParameters(x);


stimulus = data.stimulus;
time = data.time;
Rguess = 0*stimulus;

if isvector(stimulus)
	[~,Rguess]=XJWNeuronWrapper(time,stimulus,p);
else
	for i = 1:size(stimulus,2)
		[~,Rguess(:,i)]=XJWNeuronWrapper(time,stimulus(:,i),p);
	end
	clear i
end


