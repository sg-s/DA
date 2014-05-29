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
function [NLNFit, x] = FitNLNModel(data)



	x0 = [108  17  0.4  0.7   2    50   2     8000     8000   1];
	lb = [50   1    0    0.1  0    1    0     0        0      1];
	ub = [200  30   10   100  10   1e8   10     10000    10000  10];

	psoptions = psoptimset('UseParallel',true,'CompletePoll', 'on', 'Vectorized', 'off','Display','iter','MaxIter',1000,'MaxFunEvals',5000);
	x = patternsearch(@(x) NLNCostFunction(x,data,@Cost2),x0,[],[],[],[],lb,ub,psoptions);	

	% x = fminsearch(@(x) NLNCostFunction(x,data,@Cost2),x0);
	% get the final output
	NLNFit = SolveNLNModel(x,data);

function [Rguess] = SolveNLNModel(x,data)	
	% unpack data
	stimulus = data.stimulus;
	response = data.response;

	% unpack parameters of model
	tau1 = x(4); n1 = x(5); tau2 = x(6); n2 = x(7);   % bi-lobed filter block

	% pass stimulus through input non-linearity
	a = hill(x(1:3),stimulus);

	% pass output through filter
	t = 1:1000;
	K = make_bilobe_filter(tau1,n1,tau2,n2,t);
	a = filter(K,1,a);

	
	% pass filtered output through output non-linearity
	Rguess = hill(x(8:10),a(:));

function [cost] = NLNCostFunction(x,data,CostFunctionHandle)
	% unpack data
	response = data.response;

	[Rguess] = SolveNLNModel(x,data);
	% calculate the cost
	cost = CostFunctionHandle(response,Rguess);


function [r] = hill(x,xdata)
	r = x(1)*(1./(1+(x(2)./xdata).^x(3)));


function f = make_bilobe_filter(tau1,n1,tau2,n2,t)
	f1 = t.^n1.*exp(-t/tau1); % functional form in paper
	f1 = f1/tau1^(n1+1)/gamma(n1+1); % normalize appropriately


	f2 = t.^n2.*exp(-t/tau2); % functional form in paper
	f2 = f2/tau2^(n2+1)/gamma(n2+1); % normalize appropriately

	f = f1 - f2;
	f = f/max(f);