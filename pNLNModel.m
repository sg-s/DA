% pNLNModel.m
% parametric NLN model, which can by fit by FitModel2Data
% 
% created by Srinivas Gorur-Shandilya at 2:55 , 21 July 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [r] = pNLNModel(stim,p)

% all parameters of this model, for clarity
p.n1;
p.k1;  % input non-linearity

p.A;
p.tau1;
p.tau2; % filter
p.n; 
p.delay; % filter delay
p.offset;

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
lb.delay = 0;

% upper bounds
ub.A = 1;
ub.tau1 = 200;
ub.tau2 = 201;
ub.n = 3;
ub.n1 = 3;
ub.n2 = 3;
ub.delay = 100;


% input nonlinearity 
a = hill([1 p.k1 p.n1],stim); % no need for a scale parameter because it is redundant.


% choose the filter length wisely
t_max = max([p.tau1 p.tau2]);
K = filter_gamma2(1:5*t_max,p);

% is there a delay? 
p.delay = floor(p.delay);
if p.delay > 0
	K = [zeros(1,p.delay) K];
end
K = K(:);

b = filter(K,1,a);

% add an offset
b = b + p.offset;

% output nonlinearity 
r = hill([p.B p.k2 p.n2],b);
