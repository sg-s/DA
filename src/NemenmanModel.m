% NemenmanModel
% 
% 
% created by Srinivas Gorur-Shandilya at 1:49 , 16 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [R] = NemenmanModel(S,p)

% list parameters for clarity:

% input nonlinearity:
p.k; % steepness
p.x0; 

% ode:
p.d; % degradation rate 

p.lag;

% bounds
lb.k = eps;
lb.d = eps;
lb.A = eps;

ub.x0 = 0;

R = 0*S;
R(1) = 10;

S = circshift(S,round(p.lag));

for i = 2:length(S)
	dr = p.A*logistic(S(i-1),1,p.k,p.x0) - p.d*R(i-1);
	R(i) = R(i-1) + dr;
	R(i) = max([R(i) 0]);
end