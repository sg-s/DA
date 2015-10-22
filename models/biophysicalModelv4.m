% biophysicalModelv4.m
% reduced version of v3
% Euler-integration impelemtation of an abtract model implementing receptor modification
% there are two essential axes in this model: fraction of bound receptors, and fraction of methylated receptors
% 
% 
% created by Srinivas Gorur-Shandilya at 6:24 , 15 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [resp,b,m] = biophysicalModelv4(S,p)

    %     sb: 12520
    %     sm: 448
    %     kb: 0.3949
    %     km: 0.0414
    % hill_A: 24420
    % hill_k: 14.6910
    %  theta: 0.5367


% time parameters 
p.sb;
p.sm;


% ratio parameters
p.kb;
p.km;
p.theta;

% hill function output
p.hill_A;
p.hill_k;


% some bounds
lb.kb = 0;
lb.km = 0;
lb.d0 = 0;
lb.sb = 0;
lb.sm = 0;

ub.sb = 1e5;
ub.sm = 1e4;
ub.km = 10;
ub.hill_A = 1e3;

lb.theta = 0;
ub.theta = 1;



% set initial conditions
b = .1 + 0*S; 
m = .1 + 0*S;

% divide time parameters by timestep
dt = 1e-3; % hardcoded for now
p.sb = p.sb*dt;
p.sm = p.sm*dt;

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
	fx = (1-m(i-1))*S(i-1)*p.km*p.sm - p.sm*m(i-1);
	m(i) = m(i-1) + dt*fx;
	if m(i) > 1
		m(i) = 1;
	end
	if m(i) < 0
		m(i) = 0;
	end

end

% now convert bound fractioninto a firing rate using a hill function
resp = hill([p.hill_A p.hill_k 2],b);

