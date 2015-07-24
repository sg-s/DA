% pDeliverySystem.m
% NLN-based model of Delivery System
% 
% created by Srinivas Gorur-Shandilya at 2:55 , 21 July 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function [r] = pDeliverySystem(stim,p)

if nargin == 1
	% try to get the parameters from the global workspace
	global pDS
	p = pDS;
end

% all parameters of this model, for clarity
p.n1;
p.k1;  % input non-linearity

p.A;
p.tau;
p.n; 
p.offset0;
p.offset1;
p.offset2;

% output non-linearity -- depends on stimulus
p.B0;
p.B1; 
p.B2;
p.m0;
p.m1;
p.q0;
p.q1;

% lower bounds
lb.n1 = 1;
lb.k1 = eps;

lb.A = 0;
lb.tau = 1;
lb.n = 1;

lb.B0 = eps;
lb.B1 = eps;
lb.m0 = eps;
lb.m1 = eps;
lb.q0 = eps;
lb.q1 = eps;


% upper bounds
ub.n1 = 5;

ub.A = 500;
ub.tau = 300;
ub.n = 5;

ub.m0 = 3;
ub.m1 = 3;


% input nonlinearity 
a = hill([1 p.k1 p.n1],stim); % no need for a scale parameter because it is redundant.


% choose the filter length wisely
K = filter_gamma(1:5*p.tau,p);
K = K(:);

b = filter(K,1,a);

% add an offset
b = b + p.offset0 + p.offset1*mean(stim) + p.offset2*mean(stim)*mean(stim);

% output nonlinearity 
eff_B = p.B0 + mean(stim)*p.B1 + mean(stim)*mean(stim)*p.B2;
eff_M = p.m0 + mean(stim)*p.m1;
eff_Q = p.q0 + mean(stim)*p.q1;
r = hill([eff_B eff_M eff_Q],b);
