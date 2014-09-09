% SolveNLNModel2.m
% reduced version of a NLN model, with 5 parameters defining the static non-linearities
% in this model, a non-parametric filter connecting the two non-linearities is calculated on the fly. 
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function [Rguess, K, a, b] = SolveNLNModel2(p,stimulus,response)	

% pass stimulus through input non-linearity
a = hill2(p(1:2),stimulus);

% generate b(t) by inverting the output NL function:
b = ihill2(p(3:4),response);

% fit a filter from a -> b
K = FitFilter2Data(stimulus,response);

% pass a through filter
bcap = filter(K,1,a-mean(a)) + mean(b(~isinf(b)));

% pass it through a rectifier
bcap(bcap<0) = 0;

% pass filtered output through output non-linearity
Rguess = hill2(p(3:4),bcap);

