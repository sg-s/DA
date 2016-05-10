% methylationModelX.m
% simple model of receptor binding and unbinding +
% diffusible factor that decreases the rate of binding 
% diffusible factor depends on bound fraction only, not on stimulus
% 
% created by Srinivas Gorur-Shandilya at 11:19 , 26 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [b,m] = methylationModelX(S,p)

% parameters
p.r_b; % binding rate
p.r_m; % methylation rate

% ratio constants
p.k_b;
p.k_m;

% this governs the ratio of methylated to un-methylated reaction rates
p.theta;

% lower bounds
lb.r_b = 0;
lb.r_m = 0;
lb.k_b = 0;
lb.k_m = 0;
lb.theta = 0;

% upper bounds
ub.r_b = 100;
ub.r_m = 100;
ub.k_b = 100;
ub.k_m = 100;
ub.theta = 1;


b = 0*S;
m = 0*S;

for i = 2:length(S)
	M = m(i-1);
	B = b(i-1);

	fx = p.r_b*(1 - M + M*p.theta)*( (1-B*p.k_b*S(i)) - B);

	b(i) = fx + B;

	if b(i) < 0
		b(i) = 0;
	end

	if b(i) > 1
		b(i) = 1;
	end

	fx = p.r_m*(p.k_m*(1-M) - M);

	m(i) = M + fx;

	if m(i) < 0
		m(i) = 0;
	end

	if m(i) > 1
		m(i) = 1;
	end


end



