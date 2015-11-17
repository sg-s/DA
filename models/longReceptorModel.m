% longReceptorModel.m
% simple model of receptor binding and unbinding +
% diffusible factor that increases the rate of unbinding 
% 
% created by Srinivas Gorur-Shandilya at 11:19 , 26 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [R,b,d] = longReceptorModel(S,p)

% parameters
p.r_b; % binding rate
p.r_d; % diffusible factor rate

% ratio constants
p.theta_b;
p.theta_d;

p.k;

% output linear mapping 
p.A;
p.R0;

% bounds
lb.r_b = 0;
lb.r_d = 0;
lb.theta_b = 0;
lb.theta_d = 0;
ub.theta_b = 20;
ub.theta_d = 20;

lb.A = 0;
lb.k = 0;

b = 0*S; b(1) = .2;
d = 0*S; d(1) = 1e-2;

for i = 2:length(S)
	fx = (1-b(i-1))*p.r_b*S(i) - b(i-1)*p.r_b*p.theta_b*(d(i-1));
	b(i) = fx + b(i-1);

	if b(i) < 0
		b(i) = 0;
	end

	if b(i) > 1
		b(i) = 1;
	end

	fx = p.r_d*((1-b(i-1))/(b(i-1)))*exp(b(i-1)/p.k) - p.r_d*p.theta_d*d(i-1);

	d(i) = d(i-1) + fx;

	if d(i) < 0
		d(i) = 0;
	end


end

R = p.A*b + p.R0;



