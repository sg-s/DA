function [R,y,z] = DA_integrate3(S,p)
%% function [R,y,z] = DA_integrate3(S,p)
% function takes argument of stimulus and a parameter set and implements the DA
% model as described in Clark et al., 2013
% here, the model is modified as follows:
% we introduce a new parameter: m, that is a scaling parameter multiplied to the stimulus
% 								offset, a stimulus offset parameter that is added to m*S

if ~nargin
	help DA_integrate2
	return
end

% fix the stimulus
S = p.m*S + p.offset;

t = 0:1000; % filters to be this long; don't use n*tau longer than a few hundred ms in this case...
% Kz and Ky are the filters described in equations 12 and 13
Ky = generate_simple_filter(p.tau_y,p.n_y,t);
Kz = p.C*Ky + (1-p.C) * generate_simple_filter(p.tau_z,p.n_z,t);

% y and z are the stimulus convolved with the filters Ky and Kz
y = filter(Ky,1,S);
z = filter(Kz,1,S);

R = p.A*y./(1+p.B*z); % + p.r0;



function f = generate_simple_filter(tau,n,t)
f = t.^n.*exp(-t/tau); % functional form in paper
f = f/tau^(n+1)/gamma(n+1); % normalize appropriately




