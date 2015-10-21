% biophysicalModelEulerv3.m
% Euler-integration impelemtation of an abtract model implementing receptor modification
% there are two essential axes in this model: fraction of bound receptors, and fraction of methylated receptors
% 
% 
% created by Srinivas Gorur-Shandilya at 6:24 , 15 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [resp,b,m,d] = biophysicalModelv3(S,p)

    %     sb: 12520
    %     sm: 448
    %      a: 4993
    %      b: 1
    %     kb: 0.3949
    %     km: 0.0414
    %     d0: 229.0375
    % hill_A: 24420
    % hill_n: 1.4833
    % hill_k: 14.6910
    %  theta: 0.5367


% time parameters 
p.sb;
p.sm;
p.a;
p.b;

% ratio parameters
p.kb;
p.km;
p.d0;
p.theta;

% hill function output
p.hill_A;
p.hill_n;
p.hill_k;


% some bounds
lb.kb = 0;
lb.km = 0;
lb.d0 = 0;
lb.sb = 0;
lb.sm = 0;
lb.a = 0;
lb.b = 0;

ub.sb = 5e4;
ub.sm = 5e3;
ub.a = 5e3;
ub.b = 5e3;

lb.theta = 0;
ub.theta = 1;



% set initial conditions
b = 0*S; 
m = 1+0*S;
d = 1+0*S;

% divide time parameters by timestep
dt = 1e-3; % hardcoded for now
p.sb = p.sb*dt;
p.sm = p.sm*dt;
p.a = p.a*dt;
p.b = p.b*dt;

for i = 2:length(S)
	% bound fraction
	fx = (1-b(i-1))*S(i)*p.kb*p.sb*(1 + m(i-1)*(p.theta - 1) ) - p.sb*b(i-1);
	b(i) = b(i-1) + dt*fx;
	if b(i) > 1
		b(i) = 1;
	end
	if b(i) < 0
		b(i) = 0;
	end

	% methylated fraction
	fx = (1-m(i-1))*d(i-1)*p.km*p.sm - p.sm*m(i-1);
	m(i) = m(i-1) + dt*fx;
	if m(i) > 1
		m(i) = 1;
	end
	if m(i) < 0
		m(i) = 0;
	end

	% diffusible factor
	fx = p.a*b(i-1) - p.b*d(i-1) + p.d0;
	d(i) = d(i-1) + dt*fx;
	if d(i) < 0
		d(i) = 0;
	end

end

% now convert bound fractioninto a firing rate using a hill function
resp = hill([p.hill_A p.hill_k p.hill_n],b);

