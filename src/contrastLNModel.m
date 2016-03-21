% contrastLNModel
% the stimulus has two vectors:
% the first is the linear prediction
% and the second is the derivative of the stimulus
% 
% created by Srinivas Gorur-Shandilya at 1:49 , 16 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [R,Kg] = contrastLNModel(S,p)

p.n0;
p.tau;
p.K;
p.A;
p.B;
p.n;

% bounds
lb.tau = 1;
lb.K = 0;
lb.A = 1;
lb.n0 = 0;
lb.B = 0;
lb.n = 1;

ub.n = 5;
ub.tau = 300;

fp = S(:,1);
S = S(:,2);
S(S<0) = 0;

filter_length = 4*(p.n*p.tau);
t = 0:filter_length; 
Kg = generate_simple_filter(p.tau,p.n,t);

Shat = filter(Kg,1,S);

n = p.n0./(1 + p.B.*Shat);

R = real((fp.^n)./(fp.^n + p.K.^n));
R(R<0) = 0;
R = R*p.A;


function f = generate_simple_filter(tau,n,t)
f = t.^n.*exp(-t/tau); % functional form in paper
f = f/tau^(n+1)/gamma(n+1); % normalize appropriately



