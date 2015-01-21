% LN Model with gain adaptation 
% the response filter is modelled by a double gamma filter (two lobed)
% and the gain filter is modelled by a alpha filter (1 lobed)
% p.tau1 = 10;
% p.K_n = 2;
% p.tau2 = 20;
% p.K_A = 1;
% p.A = 100; % these are the parameters for the hill function at the output
% p.n = 2;
% p.Kd = 5;
% p.offset = 20;

% % gain filter parameters:
% p.tau_g = 20;
% p.beta = 1;
%
% created by Srinivas Gorur-Shandilya at 2:43 , 17 December 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [f,K,shat] = LNAModel(s,p)

% make the filters
t = 1:300;
K = filter_gamma2(p.tau1,p.K_n,p.tau2,p.K_A,t);

% filter the input
shat = filter(K,1,s-mean(s));

% add the offset
shat = shat + p.offset;

% pass through the non-linearity
x = [p.A p.Kd p.n];
f = hill(x,shat);

Kg = filter_alpha(p.tau_g,1);

f = f./(1+ p.beta*filter(Kg,1,shat-p.offset));