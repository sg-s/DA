% delivery_system.m
% part of the DA project
% detailed model of the delivery system
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function [r] = delivery_system2(f,p)

if ~nargin
	help delivery_system2
	return
end

% list parameters
p.A;
p.B;
p.C;
p.tau_y;
p.tau_z;
p.n_y;
p.n_z;
p.s0;
p.HA;
p.HB;
p.HC;


% specify bounds for FitModel2Data
lb.A = 1; lb.B = 1; lb.C = 0 ; 
lb.tau_y = 1; lb.tau_z = 1;

ub.C = 1; 
lb.HC = 0;
lb.HA = 0;
lb.HB = 0;

% extra bounds
lb.n_y = 1e-3; lb.n_z = 1e-3;
ub.n_y = 20; ub.n_z = 20;
lb.s0 = -10; ub.s0 = 10;
ub.tau_z = 400; ub.tau_y = 200;


r = DAModelv2(f,p);

r = hill([p.HA p.HB p.HC],r);


