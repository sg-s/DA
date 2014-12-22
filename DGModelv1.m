% Dynamic Gain Model v1 (LN Model-inspired).m
% see:
% https://github.com/sg-s/DA/wiki/List-of-Models
% for docs
% 
% created by Srinivas Gorur-Shandilya at 10:52 , 08 December 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function [f,shat] = DGModelv1(s,p)
switch nargin
case 0
	help DGModelv1
	return
case 1
	error('No parameters specified')
end

% make the filters
t = 1:300;
Kr = filter_gamma2(p.tau1,p.K_n,p.tau2,p.K_A,t);
Kos = filter_gamma(p.OS_tau,p.OS_n,1,t);
Kor = filter_gamma(p.OR_tau,p.OR_n,1,t);

% run through response filter
shat = filter(Kr,1,s-mean(s));

% add a offset
shat = shat + p.offset;

% compute the time-dependent Kd
stimulus_driven_Kd = p.gamma*filter(Kos,1,s);
response_driven_Kd = (1-p.gamma)*filter(Kor,1,s);
Kd = p.Kd + stimulus_driven_Kd + response_driven_Kd;

% pass through output non-linearity
f = 0*shat;
for i = 1:length(shat)
	f(i) = hill([p.A  Kd(i) p.n],shat(i));
end









