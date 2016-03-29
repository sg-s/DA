% NemenmanModel
% 
% 
% created by Srinivas Gorur-Shandilya at 1:49 , 16 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [R] = NemenmanModel(S,p)

% list parameters for clarity:
p.d; % degradation rate 
p.A;
p.lag;

% bounds
lb.d = eps;
lb.A = eps;

R = 0*S;
R(1) = 100;

% is there a lag?
S = circshift(S,round(p.lag));

% pass the stimulus through a step function
S2 = S;
ms = nanmean(S);
S2(S<ms) = 0;
S2(S>ms) = 1;
S2(S==ms) = .5;

for i = 2:length(S)
	dr = p.A*S(i-1) - p.d*R(i-1);
	R(i) = R(i-1) + dr;
	R(i) = max([R(i) 0]);
end