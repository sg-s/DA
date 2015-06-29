% iff_ifb.m
% part of the DA project
% simulates and solves using the combined IFF-IFB model as shown in pg 8 of Schulze et al. eLife 2015
% 
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function [y,u] = iff_ifb(x,p)
if ~nargin
	help iff_ifb
	return
end

% set bounds


time = 1e-3*(1:length(x));
x = x(:);
y = 0*x;

% initial conditions
y(1) = rand;
u(2) = rand;

dt = (time(2) - time(1));

for i = 2:length(time)


	% calculate derivative
	fu = p.alpha1*x(i-1) - p.alpha2*u(i-1) + p.alpha3*y(i-1);

	f1 = p.beta1*(x(i-1)/(p.beta2 + x(i-1) + p.beta3*u(i-1)));
	num = y(i-1)^p.n;
	denom = (num + p.theta^p.n);
	f2 = -p.beta4*((num)/(denom));
	f3 = -p.beta5*y(i-1);

	fy = f1+f2+f3;

	% linear extrapolation
	y(i) = y(i-1) + fy*dt;
	u(i) = u(i-1) + fu*dt;
end

