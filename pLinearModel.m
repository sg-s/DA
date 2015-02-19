% Parametric Linear Model
% the filter is modelled by a double gamma filter (two lobed)
% p.tau1 = 10;
% p.K_n = 2;
% p.tau2 = 20;
% p.K_A = 1;
% p.offset = 20;
% p.scale = 1;
% created by Srinivas Gorur-Shandilya at 2:43 , 17 December 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [f,K,shat] = pLinearModel(s,p)

% make the filters
t = 1:300;
K = filter_gamma2(p.tau1,p.K_n,p.tau2,p.K_A,t);

% filter the input
shat = filter(K,1,s-mean(s));

% add the offset
shat = shat + p.offset;

% scale
f = shat*p.scale;

f(f<0)=0;
