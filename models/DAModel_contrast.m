% DAModel_contrast.m
% modified version of the DA model that attempts to account for contrast sensation too
% 
%
% This work is licensed under a Creative Commons Attribution-ShareAlike 3.0 Unported License.
% CC-BY-SA
% Damon A. Clark, 2013


function [R,Gm,Gs] = DAModel_contrast(S,p)

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

% list parameters for legibility 
p.A;
p.B;
p.C;
p.tau_r;
p.tau_m;
p.tau_s;

lb.C = 0;
lb.A = 0;
lb.B = 0;
lb.tau_r = 0;
lb.tau_m = 0;
lb.tau_s = 0;

ub.tau_r = 1000;
ub.tau_m = 1000;
ub.tau_s = 1000;
ub.C = 2;

S = (S + p.s0); 

% see https://github.com/sg-s/DA/issues/114 for an explanation of the following
filter_length = 4*max([2*p.tau_m  2*p.tau_s 2*p.tau_r]);
if filter_length < length(S)/10
else
	filter_length = length(S)/10; % ridiculously long filters
end
t = 0:filter_length; 

% there are 3 filters: Kr, the response filter, Km, the mean gain control filter, and Ks, the variance gain control filter
Kr = generate_simple_filter(p.tau_r,2,t);
Km = generate_simple_filter(p.tau_m,2,t);
Ks = generate_simple_filter(p.tau_s,2,t);


% generate response by convolving stimulus with response filter
R = filter(Kr,1,S);

% now calculate the gain
Gm = filter(Km,1,S); 
Gs = filter(Ks,1,S); 
Gs = S - Gs;
Gs(Gs<0) = 0;
G = 1+ p.B*(Gm + p.C*Gs);

R = p.A*R./G;




function f = generate_simple_filter(tau,n,t)
f = t.^n.*exp(-t/tau); % functional form in paper
f = f/tau^(n+1)/gamma(n+1); % normalize appropriately




