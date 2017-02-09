% contrastGainBeta.m
% 
% created by Srinivas Gorur-Shandilya at 2:12 , 10 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


function R = contrastGainBeta(S,p)

% only one parameter
p.A;
p.B;

lb.B = 0;
lb.A = 1;
ub.B = 1e4;
ub.A  = 1e4;

fp = S(:,2);
S = S(:,1);

R = p.A*fp./(1+p.B*S);