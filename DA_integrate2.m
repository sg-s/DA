function [R,y,z] = DA_integrate2(S,p)

%% function takes argument of stimulus and a parameter set and implements the DA
%% model as described in Clark et al., 2013
%% the parameter set is as described in DA_model_script.m

%% This work is licensed under a Creative Commons Attribution-ShareAlike 3.0 Unported License.
%% CC-BY-SA
%% Damon A. Clark, 2013
% 
% modified by Srinivas Gorur-Shandilya
% here we throw away tau_r as we never use it and instead introduce a new parameter, r0, which we add on top to the response at the end. 


t = [0:3000]; % filters to be this long; don't use n*tau longer than a few hundred ms in this case...
% Kz and Ky are the filters described in equations 12 and 13
Ky = generate_simple_filter(p.tau_y,p.n_y,t);
Kz = p.C*Ky + (1-p.C) * generate_simple_filter(p.tau_z,p.n_z,t);

% y and z are the stimulus convolved with the filters Ky and Kz
y = filter(Ky,1,S);
z = filter(Kz,1,S);

R = p.A*y./(1+p.B*z) + p.r0;



function f = generate_simple_filter(tau,n,t)
f = t.^n.*exp(-t/tau); % functional form in paper
f = f/tau^(n+1)/gamma(n+1); % normalize appropriately



