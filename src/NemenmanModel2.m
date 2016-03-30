% NemenmanModel2.m
% here we extend the Nemenman contrast-invariant model with a filter
% 
% created by Srinivas Gorur-Shandilya at 1:49 , 16 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [R] = NemenmanModel2(S,p)

% list parameters for clarity:
p.d; % degradation rate 
p.B;
% filter parameters 
p.n;
p.tau1;
p.tau2;
p.A;

% bounds
lb.d = eps;
lb.A = eps;
lb.B = eps;
lb.tau1 = 1;
lb.tau2 = 2;
lb.n = 2;

ub.n = 2;
ub.tau1 = 300;
ub.tau2 = 500;

R = 0*S;
R(1) = 10; % initial condition, might have to dick around with this to have nice solutions. 

% pass the stimulus through a step function
S2 = S;
ms = nanmean(S);
S2(S<ms) = 0;
S2(S>ms) = 1;
S2(S==ms) = .5;
S = S2; clear S2

for i = 2:length(S)
	dr = p.B*S(i-1) - p.d*R(i-1);
	R(i) = R(i-1) + dr;
	R(i) = max([R(i) 0]);
end

% now pass it through a filter
K = filter_gamma2([],p);

R = filter(K,1,R);