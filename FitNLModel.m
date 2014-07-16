% FitNLModel.m
% fits a Nonlinear-Linear Model to data
% usage: [NLFit p] = FitNLModel(data)
% where data is a structure containing 2 vectors:
% stimulus
% response
%
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function [NLFit, x] = FitNLModel(data,x0)

switch nargin 
	case 0
		help FitNLModel
		return
	case 1
		x0 = [40   5    1200];
		lb = [5    -20    800  ];
		ub = [200   30   2000 ];

	case 2
		lb = x0/10;
		ub = x0*10;
end


psoptions = psoptimset('UseParallel',true,'CompletePoll', 'on', 'Vectorized', 'off','Display','iter','MaxIter',100,'MaxFunEvals',10000);
x = patternsearch(@(x) NLCostFunction(x,data),x0,[],[],[],[],lb,ub,psoptions);	

% get the final output
NLFit = SolveNLModel(x,data.stimulus,data.response);

% debug
figure, hold on
plot(data.response), hold on
plot(NLFit,'r')