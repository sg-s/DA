% DAModel.m
% modified to have a bilobed Kz filter. 
%
% This work is licensed under a Creative Commons Attribution-ShareAlike 3.0 Unported License.
% CC-BY-SA
% Damon A. Clark, 2013
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [R,y,z,Ky,Kz] = DAModel_bilobedKz(S,p)

switch nargin
case 0
	help DA_Modelv2
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

% specify bounds for FitModel2Data
lb.A = 1; lb.B = 0; lb.C = 0 ; 
lb.tau_y = 1; lb.tau_z1 = 1; lb.tau_z2 = 1;

ub.C = 1; 

% extra bounds
lb.n_y = 2; lb.n_z = 2;
ub.n_y = 2; ub.n_z = 2;
lb.s0 = -5; ub.s0 = 1;
ub.tau_z1 = 400; ub.tau_y = 200;
ub.tau_z2 = 400; 



S = (S + p.s0); 

% see https://github.com/sg-s/DA/issues/114 for an explanation of the following
filter_length = 4*max([p.n_z*p.tau_z2  p.n_y*p.tau_y]);
if filter_length < length(S)/10
else
	filter_length = length(S)/10; % ridiculously long filters
end
t = 0:filter_length; 

% Kz and Ky are the filters described in equations 12 and 13
q.A = p.Az;
q.n = p.n_z;
q.tau1 = p.tau_z1;
q.tau2 = p.tau_z2;
Ky = generate_simple_filter(p.tau_y,p.n_y,t);
Kz = p.C*Ky + (1-p.C) * filter_gamma2(t,q);

% y and z are the stimulus convolved with the filters Ky and Kz
y = filter(Ky,1,S);
z = filter(Kz,1,S);

R = (p.A*y./(1+p.B*z));




function f = generate_simple_filter(tau,n,t)
f = t.^n.*exp(-t/tau); % functional form in paper
f = f/tau^(n+1)/gamma(n+1); % normalize appropriately




