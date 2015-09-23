function [R,y,z] = DA_integrate3(S,p,LNpred)
%% function [R,y,z] = DA_integrate2(S,p,LNpred)
% function takes argument of stimulus and a parameter set and implements the DA
% model as described in Clark et al., 2013
% modified as follows: it directly modulates the LN prediction, which it accepts as a vector
% the parameter set is as described in DA_model_script.m
%
% This work is licensed under a Creative Commons Attribution-ShareAlike 3.0 Unported License.
% CC-BY-SA
% Damon A. Clark, 2013
% 
% modified by Srinivas Gorur-Shandilya
% here we throw away tau_r as we never use it and instead introduce a new parameter, s0, which is subtracted from the stimulus. 

if ~nargin
	help DA_integrate3
	return
end

S = (S - p.s0); 

t = [0:1000]; % filters to be this long; don't use n*tau longer than a few hundred ms in this case...
% Kz and Ky are the filters described in equations 12 and 13
Ky = generate_simple_filter(p.tau_y,p.n_y,t);
Kz = p.C*Ky + (1-p.C) * generate_simple_filter(p.tau_z,p.n_z,t);

% y and z are the stimulus convolved with the filters Ky and Kz
y = filter(Ky,1,S);
z = filter(Kz,1,S);

% R = (p.A*y./(1+p.B*z));
R = ((LNpred)./(1+p.B*z));

% % pass through rectifier
% R(R<0) = 0;


function f = generate_simple_filter(tau,n,t)
f = t.^n.*exp(-t/tau); % functional form in paper
f = f/tau^(n+1)/gamma(n+1); % normalize appropriately




