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
function [NLNFit, x, K] = FitNLNModel(data,x0)

switch nargin 
	case 0
		help FitNLNModel
		return
	case 1
		x0 = [4    1         1    6   ];
		lb = [1    1         1    1   ];
		ub = [200  30        10   10  ];

	case 2
		lb = x0/10;
		ub = x0*10;
end


psoptions = psoptimset('UseParallel',true,'CompletePoll', 'on', 'Vectorized', 'off','Display','iter','MaxIter',200,'MaxFunEvals',10000);
x = patternsearch(@(x) NLNCostFunction(x,data),x0,[],[],[],[],lb,ub,psoptions);	

% get the final output
[NLNFit, K] = SolveNLNModel2(x,data.stimulus,data.response);


% % debug
% figure, hold on
% plot(data.response,'b')
% plot(NLNFit,'r')