% simpleReceptorModelX.m
% simple model of receptor binding and unbinding +
% diffusible factor that decreases the rate of binding 
% diffusible factor depends on bound fraction only, not on stimulus
% 
% created by Srinivas Gorur-Shandilya at 11:19 , 26 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [R,b,d] = simpleReceptorModelX(S,p)

% parameters
p.r_b; % binding rate
p.r_d; % diffusible factor rate

% ratio constants
p.theta_b;
p.theta_d;

% output nonlinearity 
p.hill_A;
p.hill_K;

% bounds
lb.r_b = 0;
lb.r_d = 0;
lb.theta_b = 0;
lb.theta_d = 0;
ub.theta_b = 10;
ub.theta_d = 10;

lb.hill_A = 0;
lb.hill_K = 0;
ub.hill_K = 1;


b = 0*S;
d = 0*S;

for i = 2:length(S)
	fx = (1-b(i-1))*p.r_b*S(i)/(1+d(i-1)) - b(i-1)*p.r_b*p.theta_b;
	b(i) = fx + b(i-1);

	if b(i) < 0
		b(i) = 0;
	end

	if b(i) > 1
		b(i) = 1;
	end

	fx = p.r_d*b(i-1) - p.r_d*p.theta_d*d(i-1);

	d(i) = d(i-1) + fx;

	if d(i) < 0
		d(i) = 0;
	end


end


R = p.hill_A*(b.^2)./(b.^2 + p.hill_K^2);




