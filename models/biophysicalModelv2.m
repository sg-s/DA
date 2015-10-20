% biophysicalModelEulerv2.m
% Euler-integration impelemtation of a reduced version of the Nagel-Wilson Model
% this version differs from the older one in that we use a parametric LN model to transform channel opening to firing rates. 
% 
% created by Srinivas Gorur-Shandilya at 6:24 , 15 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [resp,receptors,channels,diff_factor] = biophysicalModelv2(S,p)

% parameters concerning the receptor model:
p.kb;
p.sb;
p.ko;
p.kc;
p.a;
p.b;

% parameters about the output nonlinearity 
p.hill_A;
p.hill_n;
p.hill_k;


% some basic bounds
lb.kb = 0;
lb.sb = 0;
lb.ko = 0;
lb.kc = 0;
lb.a = 0;
lb.b = 0;

lb.hill_A = 0;
lb.hill_n = 1;
lb.hill_k = -1;

ub.kb = 1e3;
ub.sb = 1e3;
ub.ko = 1e3;
ub.kc = 1e3;
ub.a = 1e3;
ub.b = 1e3;

ub.n = 10;


% 3 dimensional matrix, make some placeholders
x = NaN(3,length(S)); % 6 dimensions as listed in the paper + 1 for LFP voltage

% define initial conditions 
x(1,1) = .5; % half of all receptors bound
x(2,1) = 0.5; % half of all channels open
x(3,1) = 0; % no diffusable factor 


dt = 1e-3; % hard coded for now

% divide time dependant parameters by the time step
p.kb = p.kb*dt;
p.a = p.a*dt;
p.b = p.b*dt;


for i = 2:length(S)
	O = S(i-1); % odor stimulus

	% get the receptor activities right
	fx = (1-x(1,i-1))*O*p.kb*p.sb - x(1,i-1)*p.sb;
	x(1,i) = x(1,i-1) + dt*fx;

	if x(1,i) < 0
		x(1,i) = 0;
	end
	if x(1,i) > 1
		x(1,i) = 1;
	end

	% now the channel opening 
	fx = p.ko*(1-x(2,i-1))*x(1,i-1) - x(2,i-1)*x(3,i-1)*p.kc; % modified from Nagel and Wilson
	x(2,i) = x(2,i-1) + dt*fx;

	if x(2,i) < 0
		x(2,i) = 0;
	end
	if x(2,i) > 1
		x(2,i) = 1;
	end


	% now the diffusable factor
	fx = p.a*x(2,i-1) - p.b*x(3,i-1);  
	x(3,i) = x(3,i-1) + dt*fx;
	if x(3,i) < 0
		x(3,i) = 0;
	end

end	

diff_factor = x(3,:);
channels = x(2,:);
receptors = x(1,:);

% now convert channel opening probability into a firing rate using a hill function
resp = hill([p.hill_A p.hill_k p.hill_n],channels);

