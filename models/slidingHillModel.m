% slidingHillModel
% a model that modulates the K parameter of a Hill function in a manner
% dependent on the stimulus
% 
% created by Srinivas Gorur-Shandilya at 1:49 , 16 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [R,Kg,K] = slidingHillModel(S,p)

% for stimulus
p.x_offset;

% for Hill function
p.A;
p.n;

% for changing K:
p.K0;
p.m;
p.tau;
p.B;

% bounds
lb.tau = 1;
lb.n = 1;
lb.m = 1;
lb.B = eps;
lb.A = 10;
lb.K0 = 0;

ub.tau = 400;
ub.n = 5;
ub.m = 10;
ub.B = 1000;


% apply offset
S = S + p.x_offset;

filter_length = min([4*(p.m*p.tau) 20]);
t = 0:filter_length; 
Kg = generate_simple_filter(p.tau,p.m,t);

Shat = filter(Kg,1,S);
Shat(Shat < 0) = 0;

K = p.K0./(1 + p.B*Shat);

S(S < 0) = 0;
K(K < 0) = 0;

R = (S.^p.n)./(S.^p.n + K.^p.n);
R = R*p.A;


function f = generate_simple_filter(tau,n,t)
f = t.^n.*exp(-t/tau); % functional form in paper
f = f/tau^(n+1)/gamma(n+1); % normalize appropriately



