% FitDAModelToData.m
% fits the DA model to data
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
%
% data is a structure containing:
% .stimulus -- the stimulus
% .response -- the response of the ORN (or whatever sensor)
% .time -- a time vector
% 
% stimulus or response can be matrices, where each row represents a different experiment or a different trial. 
function [p, Rguess] = FitDAModelToData(data,x0,lb,ub,IgnoreInitial)

switch nargin 
	case 0
		help FitDAModelToData
		return
	case 1
		x0 = [400   20   0.1 0.85  2   75    2   -9];
		lb = [200   0    0   0     2   1     2   -20];
		ub = [60000 900  1   10    2   100   2   10];
		IgnoreInitial = 300;
	case 2
		lb = x0/2;
		ub = x0*10;
		IgnoreInitial = 300;
	case 3
		ub = x0*10;
		IgnoreInitial = 300;
	case 4
		IgnoreInitial = 300;
end



psoptions = psoptimset('UseParallel',true,'CompletePoll', 'on', 'Vectorized', 'off','Display','iter','MaxIter',300,'MaxFunEvals',10000);
x = patternsearch(@(x) DA_cost_function(x,data,@Cost2,IgnoreInitial),x0,[],[],[],[],lb,ub,psoptions);
p = ValidateDAParameters2(x);


stimulus = data.stimulus;
Rguess = 0*stimulus;

if isvector(stimulus)
	Rguess = DA_integrate2(stimulus,p);
else
	for i = 1:size(stimulus,2)
		Rguess(:,i) = DA_integrate2(stimulus(:,i),p);
	end
	clear i
end

Rguess = abs(Rguess);