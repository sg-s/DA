function [R,y,z] = DA_integrate2(S,p)
%% function [R,y,z] = DA_integrate2(S,p)
% function takes argument of stimulus and a parameter set and implements the DA
% model as described in Clark et al., 2013
% the parameter set is as described in DA_model_script.m
%
% This work is licensed under a Creative Commons Attribution-ShareAlike 3.0 Unported License.
% CC-BY-SA
% Damon A. Clark, 2013
% 
% modified by Srinivas Gorur-Shandilya
% here we throw away tau_r as we never use it and instead introduce a new parameter, s0, which is subtracted from the stimulus. 
switch nargin
case 0
	help DA_integrate2
	return
case 1
	error('Not enough input arguments')
case 2
	if ~isvector(S)
		error('First argument should be a vector')
	end
	if ~isstruct(p)
		error('Second argument should be a structure')
	end
end



S = (S + p.s0); 

t = 0:300; % filters to be this long; don't use n*tau longer than a few hundred ms in this case...
% Kz and Ky are the filters described in equations 12 and 13
Ky = generate_simple_filter(p.tau_y,p.n_y,t);
Kz = p.C*Ky + (1-p.C) * generate_simple_filter(p.tau_z,p.n_z,t);

% y and z are the stimulus convolved with the filters Ky and Kz
y = filter(Ky,1,S);
z = filter(Kz,1,S);

R = (p.A*y./(1+p.B*z));

% % pass through rectifier
% R(R<0) = 0;


function f = generate_simple_filter(tau,n,t)
f = t.^n.*exp(-t/tau); % functional form in paper
f = f/tau^(n+1)/gamma(n+1); % normalize appropriately




