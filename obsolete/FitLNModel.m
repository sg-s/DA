% FitLNModel.m
% fits a Linear-Nonlinear Model to data
% usage: [LNFit, x] = FitLNModel(data)
% where data is a structure containing 2 vectors:
% stimulus
% response
%
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function [LNFit, x] = FitLNModel(data,x0)

switch nargin 
	case 0
		help FitLNModel
		return
	case 1
		x0 = [ 1  2];
		lb = [0 0 ];
		ub = [10 20];

	case 2
		lb = x0/10;
		ub = x0*10;
end


psoptions = psoptimset('UseParallel',true,'CompletePoll', 'on', 'Vectorized', 'off','Display','iter','MaxIter',100,'MaxFunEvals',10000);
x = patternsearch(@(x) LNCostFunction(x,data),x0,[],[],[],[],lb,ub,psoptions);	

% get the final output
[LNFit,K] = SolveLNModel(x,data.stimulus,data.response);

% debug
figure, hold on
plot(data.response), hold on
plot(LNFit,'r')

% plot the filter
figure, plot(K)

% plot the nonlinearity 
figure, hold on
a = ihill2(x,data.response);
plot(a,data.response,'k.')



disp(rsquare(LNFit,data.response))