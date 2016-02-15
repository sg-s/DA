% DAModelv2.m
% modified version of DA model, skips integration
% [R,y,z,Ky,Kz] = DAModelv2(S,p)
% function takes argument of stimulus and a parameter set and implements the DA
% model as described in Clark et al., 2013
% the parameter set is as described in DA_model_script.m
%
%
% This work is licensed under a Creative Commons Attribution-ShareAlike 3.0 Unported License.
% CC-BY-SA
% Damon A. Clark, 2013


function [R,y,z,Ky,Kz] = DAModelv2(S,p)

switch nargin
case 0
	help DA_Modelv2
	return
case 1
	error('Not enough input arguments')
case 2
	assert(isvector(S),'First argument should be a vector')
	assert(isstruct(p),'Second argument should be a structure')
end

% specify bounds for fitModel2Data
lb.A = 1; lb.B = 1; lb.C = 0 ; 
lb.tau_y = 1; lb.tau_z = 1;

ub.C = 1; 

% extra bounds
lb.n_y = 1; lb.n_z = 1;
ub.n_y = 10; ub.n_z = 5;
lb.s0 = -5; ub.s0 = 1;
ub.tau_z = 300; ub.tau_y = 200;



S = (S + p.s0); 

% see https://github.com/sg-s/DA/issues/114 for an explanation of the following
filter_length = 4*max([p.n_z*p.tau_z  p.n_y*p.tau_y]);
if filter_length < length(S)/10
else
	filter_length = length(S)/10; % ridiculously long filters
end
t = 0:filter_length; 

% Kz and Ky are the filters described in equations 12 and 13
Ky = generate_simple_filter(p.tau_y,p.n_y,t);
Kz = p.C*Ky + (1-p.C) * generate_simple_filter(p.tau_z,p.n_z,t);

% y and z are the stimulus convolved with the filters Ky and Kz
y = filter(Ky,1,S);
z = filter(Kz,1,S);

R = (p.A*y./(1+p.B*z));




function f = generate_simple_filter(tau,n,t)
f = t.^n.*exp(-t/tau); % functional form in paper
f = f/tau^(n+1)/gamma(n+1); % normalize appropriately




