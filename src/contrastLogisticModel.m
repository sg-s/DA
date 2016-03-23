% contrastLogisticModel
% the stimulus has two vectors:
% the first is the linear prediction
% and the second is the derivative of the stimulus
% 
% created by Srinivas Gorur-Shandilya at 1:49 , 16 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [R,Kg] = contrastLogisticModel(S,p)

% for logistic function
p.x0;
p.A;

% for changing k:
p.k0;
p.n;
p.tau;
p.B;

% bounds
lb.tau = 10;
lb.n = 2;
lb.B = eps;
lb.A = 1;
lb.k0 = 10;
lb.x0 = -3;

ub.tau = 400;
ub.n = 2;
ub.B = 1000;

Sd = S;
Sd(Sd<0) = 0;

filter_length = 4*(p.n*p.tau);
t = 0:filter_length; 
Kg = generate_simple_filter(p.tau,p.n,t);

Shat = filter(Kg,1,Sd);

k = p.k0./(1 + p.B*Shat);

R = logistic(S,p.A,k,p.x0);


function f = generate_simple_filter(tau,n,t)
f = t.^n.*exp(-t/tau); % functional form in paper
f = f/tau^(n+1)/gamma(n+1); % normalize appropriately



