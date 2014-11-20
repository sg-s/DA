% FitFractionalDModel.m
% fits a fractional derivative model to stimulus-response data
% the idea is simple: the response is just the fractional derivative, passed through some unknown output nonlinearity. 
% usage:
% [alpha,d,fp] =  FitFractionalDModel(stimulus, response,nsteps)
% where stimulus and response are vectors
% alpha is the degree of fractional differentation, from 0 to 1
% and d is the delay of the stimulus and the response, as determined by this script
% and NLparam is a vector containing parameters of the hill function that is fit
%  
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [alpha,d,fp] =  FitFractionalDModel(stimulus, response,nsteps)
switch nargin 
case 0
	help FitFractionalDModel
	return
case 1
	error('Not enough input arguments')
case {2,3}
	if ~isvector(stimulus) || ~isvector(response)
		error('Arguments should be vectors')
	end
	stimulus = stimulus(:);
	response = response(:);
	mr = mean(response);
	response = response/mr;
	
end
if nargin < 3
	nsteps  = 30;
end

x0 = .15;
psoptions = psoptimset('UseParallel',true, 'Vectorized', 'off','Cache','on','CompletePoll','on','Display','iter','MaxIter',nsteps,'MaxFunEvals',2000);

lb = [0];% 0 0 0];
ub = [1];% Inf Inf Inf];
x = patternsearch(@(x) FDCostFunction(x,stimulus,response),x0,[],[],[],[],lb,ub,psoptions);
alpha = x;
fp = fgl_deriv(x,stimulus,1);
fp(fp<0)= 0;
fp = fp/mean(fp);
d = finddelay(fp,response);
if d
	fp = [NaN(d,1); fp];
	fp = fp(1:length(response));
end
fp = fp*mr;


	function c= FDCostFunction(x,stimulus,response)
		fp = fgl_deriv(x,stimulus,1);
		fp(fp<0)= 0;
		fp = fp/mean(fp);
		d = finddelay(fp,response);
		if d
			c=Cost2(response(d:end),fp(1:end-d+1));
		else
			c = Cost2(response,fp);
		end

	end
end