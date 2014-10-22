% FitFractionalDModel.m
% fits a fractional derivative model to stimulus-response data
% the idea is simple: the response is just the fractional derivative, passed through some unknown output nonlinearity. 
% usage:
% [alpha, NLparam] =  FitFractionalDModel(stimulus, response); 
% where stimulus and response are vectors
% alpha is the degree of fractional differentation, from 0 to 1
% and NLparam is a vector containing parameters of the hill function that is fit
%  
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [alpha,d,fp] =  FitFractionalDModel(stimulus, response)
switch nargin 
case 0
	help FitFractionalDModel
	return
case 1
	error('Not enough input arguments')
case 2
	if ~isvector(stimulus) || ~isvector(response)
		error('Arguments should be vectors')
	end
	stimulus = stimulus(:);
	response = response(:);
	response = response/mean(response);
end


x0 = .35;
nsteps = 10;
psoptions = psoptimset('UseParallel',true, 'Vectorized', 'off','Cache','on','CompletePoll','on','Display','iter','MaxIter',nsteps,'MaxFunEvals',2000);

lb = [0];% 0 0 0];
ub = [1];% Inf Inf Inf];
x = patternsearch(@(x) FDCostFunction(x,stimulus,response),x0,[],[],[],[],lb,ub,psoptions);
alpha = x;
fp = fgl_deriv(x,stimulus,1);
fp(fp<0)= 0;
fp = fp/mean(fp);
d = finddelay(fp,response);
fp = [NaN(d,1); fp];
fp = fp(1:length(response));



	function c= FDCostFunction(x,stimulus,response)
		fp = fgl_deriv(x,stimulus,1);
		fp(fp<0)= 0;
		fp = fp/mean(fp);
		d = finddelay(fp,response);
		c=Cost2(response(d:end),fp(1:end-d+1));
	end
end