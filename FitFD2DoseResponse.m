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

function [alpha,d,fp] =  FitFD2DoseResponse(stimulus, response,a,z)
switch nargin 
case 0
	help FitFractionalDModel
	return
case {1,2,3}
	error('Not enough input arguments')
case 4
	if ~isvector(stimulus) || ~isvector(response)
		error('Arguments should be vectors')
	end
	stimulus = (stimulus(:));
	stimulus(1) = stimulus(2);
	response = response(:);

	
end

max_delay  =35;
nsteps = 20;
x0 = .35;
psoptions = psoptimset('UseParallel',false, 'Vectorized', 'off','Cache','on','CompletePoll','on','Display','none','MaxIter',nsteps,'MaxFunEvals',2000);

lb = [0];% 0 0 0];
ub = [1];% Inf Inf Inf];
x = patternsearch(@(x) FDCostFunction(x,stimulus,response),x0,[],[],[],[],lb,ub,psoptions);
alpha = x;
fp = fgl_deriv(x,stimulus,1);
m = max(response);
mprime = max(fp(a:z));
bprime = mean(fp(end-1000:end));
b = mean(response(end-1000:end));
fp = (fp - bprime)/(mprime - bprime);
fp  =fp*(m-b) + b;
fp(fp<0)=0;
fp(fp>1e3)=0;
d = finddelay(stimulus,response);
if d
	if d < max_delay 
		fp = [NaN(d,1); fp];
		fp = fp(1:length(response));
		fp = [NaN(d,1); fp];
		fp = fp(1:length(response));
	end
end

if ~nargout
	figure, hold on
	plot(response)
	hold on
	plot(fp)
end


	function c= FDCostFunction(x,stimulus,response)
		fp = fgl_deriv(x,stimulus,1);
		m = max(response);
		mprime = max(fp(a:z));
		bprime = mean(fp(end-1000:end));
		b = mean(response(end-1000:end));
		fp = (fp - bprime)/(mprime - bprime);
		fp  =fp*(m-b) + b;
		fp(fp<0)=0;
		fp(fp>1e3)=0;

		% filter 
		fp(isnan(fp)) = 0;
		fp = filtfilt(ones(1,10)/10,1,fp);

		d = finddelay(stimulus,response);
		
		if d
			if d < max_delay 
				fp = [NaN(d,1); fp];
				fp = fp(1:length(response));
				c=Cost2(response(a:z),fp(a:z));
			else
				c = Cost2(response,fp);
				% c=1-rsquare(response(a:z),fp(a:z));
			end
		else
			c = Cost2(response,fp);
			% c=1-rsquare(response(a:z),fp(a:z));
		end
		if isnan(c)
			keyboard
		end
	end
end