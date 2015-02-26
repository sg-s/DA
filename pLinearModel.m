% Parametric Linear Model
% the filter is modelled by a double gamma filter (two lobed)
% created by Srinivas Gorur-Shandilya at 2:43 , 17 December 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [f,K,shat] = pLinearModel(s,p)

% set bounds
lb.n = 0; 
ub.n = 6;
lb.tau1 = 1;
ub.tau1 = 100;
lb.tau2 = 1;
ub.tau2 = 200;
lb.A = 0;
ub.A = 2;


% show parameters for readability
p.n;
p.A;
p.tau1;
p.tau2;
p.scale;
p.offset;

% make the filters
K = filter_gamma2(p,1:300);

% filter the input
shat = filter(K,1,s-mean(s));

% add the offset
shat = shat + p.offset;

% scale
f = shat*p.scale;
