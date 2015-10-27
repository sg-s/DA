% pLinearModel.m
% parametric linear model, that accounts for some trivial offsets
% uses filter_gamma2, so you can use bounds to constrain it to a single lobed filter if needed
% 
% created by Srinivas Gorur-Shandilya at 9:57 , 26 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [R,K] = pLinearModel(S,p)

S = S(:);

% parameters
p.tau;
p.A;
p.n;

ub.A = 0;

lb.n = 1;
ub.n = 10;

lb.tau = 1;
ub.tau = 3e2;

t = 1:1000;
K = filter_gamma(t,p);

R = filter(K,1,S-mean(S));
