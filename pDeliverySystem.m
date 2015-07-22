% pDeliverySystem.m
% parametric NLN model, which can by fit by FitModel2Data
% 
% created by Srinivas Gorur-Shandilya at 2:55 , 21 July 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [r] = pDeliverySystem(stim,p)

% all parameters of this model, for clarity
p.n1;
p.k1;  % input non-linearity

p.A;
p.tau1;
p.tau2; % filter
p.n; 

p.B; 
p.n2;
p.k2; 

% lower bounds
lb.n1 = 1;
lb.k1 = eps;
lb.tau1 = 1;
lb.tau2 = 1;
lb.B = eps;
lb.n2 = 1;
lb.k2 = eps;
lb.n = 1;
lb.A = 0;

% upper bounds
ub.A = 1;
ub.tau1 = 200;
ub.tau2 = 201;
ub.n = 3;
ub.n1 = 3;
ub.n2 = 3;


% filter input
% choose the filter length wisely
t_max = max([p.tau1 p.tau2]);
K = filter_gamma2(1:5*t_max,p);
b = filter(K,1,stim);
b(b<0) = 0; % have to throw away -ve values as output nonlinearity only accepts positive values

% output nonlinearity 
r = hill([p.B mean(stim)*p.k2 mean(stim)*p.n2],b);
