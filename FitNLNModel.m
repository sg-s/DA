% FitNLNModel.m
% fits a Nonlinear-Linear-Nonlinear Model to data
% usage: [NLNFit p] = FitNLNModel(data)
% where data is a structure containing 2 vectors:
% stimulus
% response
%
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function [NLNFit, x] = FitNLNModel(data,x0)

switch nargin 
	case 0
		help FitNLNModel
		return
	case 1
		x0 = [60   1   0.4   1.5   1    191   0.1   8288     4000   1];
		lb = [50   1    0    0.1   0    1     0     0        0      0];
		ub = [200  30   10   100   10   1e8   10    50000    10000  10];

	case 2
		lb = x0/10;
		ub = x0*10;
end


psoptions = psoptimset('UseParallel',true,'CompletePoll', 'on', 'Vectorized', 'off','Display','iter','MaxIter',500,'MaxFunEvals',10000);
x = patternsearch(@(x) NLNCostFunction(x,data,@Cost2),x0,[],[],[],[],lb,ub,psoptions);	

% get the final output
NLNFit = SolveNLNModel(x,data.stimulus);
