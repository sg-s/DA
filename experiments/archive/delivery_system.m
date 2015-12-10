% delivery_system.m
% part of the DA project
% detailed model of the delivery system
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function [y,x] = delivery_system(f,p)
F = 20;
if ~nargin
	help delivery_system
	return
end

% set bounds
lb.k1 = eps;
lb.k0 = eps;
lb.k2 = eps;
lb.n = 2;
ub.n = 10;
lb.f0 = eps;
lb.k3 = eps;
lb.k5 = eps;

time = 1e-2*(1:length(f));
f = f(:);
y = 0*f;
x = 0*f;
z = 0*y;

% initial conditions
y(1) = rand;
x(2) = rand;

dt = (time(2) - time(1));

for i = 2:length(time)


	% calculate derivative
	dx = p.k1 - p.k0*x(i-1) + p.k2*hill([1 p.f0 p.n],f(i-1)) - p.k3*f(i-1)*x(i-1);

	dy = p.k3*f(i-1)*x(i-1) - p.k4*(f(i-1) + F)*y(i-1) -  p.k5*y(i-1) + p.k6*z(i-1);

	dz = - p.k6*z(i-1) + p.k5*y(i-1);

	% linear extrapolation
	y(i) = y(i-1) + dy*dt;
	x(i) = x(i-1) + dx*dt;
	z(i) = z(i-1) + dz*dt;

end

