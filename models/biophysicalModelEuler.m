% biophysicalModelEuler.m
% Euler-integration impelemtation of a reduced version of the Nagel-Wilson Model
% 
% created by Srinivas Gorur-Shandilya at 6:24 , 15 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [R,x] = biophysicalModelEuler(S,p)

% list all parameters for legibility 
p.kb;
p.sb;
p.ko;
p.kc;
p.a;
p.b;

p.stim_scale;

% bounds
lb.kb = 9000;
lb.sb = 900;
ub.kb = 11000;
ub.sb = 1100;


lb.ko = eps;
lb.kc = eps;
lb.a = eps;
lb.b = eps;

p.g_OR = eps;
p.g_leak = eps;


ub.V_rest = -50;
lb.V_OR = 10;


% 7 dimensional matrix, make some placeholders
x = NaN(4,length(S)); % 6 dimensions as listed in the paper + 1 for LFP voltage

% define initial conditions 
x(1,1) = .5; % half of all receptors bound
x(2,1) = 0.5; % half of all channels open
x(3,1) = 1; % some diffusable factor 
x(4,1) = -10; % mV??

dt = 1e-3; % hard coded for now

% divide time dependant parameters by the time step
p.ka = p.ka*dt;
p.sa = p.sa*dt;
p.kb = p.kb*dt;
p.sb = p.sb*dt;
p.a = p.a*dt;
p.b = p.b*dt;

S  = S*p.stim_scale;


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
	fx = (x(1,i-1))*p.ko*(1/x(3,i-1))*(1-x(2,i-1)) - p.kc*x(2,i-1); % according to Nagel and Wilson
	x(2,i) = x(2,i-1) + dt*fx;

	% now the diffusable factor
	fx = p.a*x(2,i-1) - p.b*x(3,i-1);  
	x(3,i) = x(3,i-1) + dt*fx;

	% now we convert channel opening to voltage using something I made up because this is where the model description ends in the paper
	fx = -(x(2,i-1)*p.g_OR*(x(4,i-1) -  p.V_OR)) - (p.g_leak*(x(4,i-1) - p.V_rest))
	x(4,i) =x(4,i-1) + dt*fx;

end	

% now convert this fake LFP into a fake firing rate
t = 1:50; % in ms
K = normpdf(t,10,3).*(-8e-3*t+.0812);
K = K/max(K);
K = K*20;

R =  (filter(K,1,x(4,:)));

R = R + p.f0;


