% biophysicalModelEuler.m
% Euler-integration impelemtation of a reduced version of the Nagel-Wilson Model
% 
% created by Srinivas Gorur-Shandilya at 6:24 , 15 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [resp,receptors,channels,diff_factor,voltage] = biophysicalModelEuler(S,p)

% list all parameters for legibility 
p.kb;
p.sb;
p.ko;
p.kc;
p.a;
p.b;
p.stim_scale;

% bounds
lb.kb = 1;
lb.sb = 0;
lb.ko = 0;
lb.kc = 0;
lb.a = 0;
lb.b = 0;

lb.g_OR = 0;
lb.g_leak = 0;

lb.V_rest = -100;
ub.V_rest = -50;
lb.V_OR = 0;
ub.V_OR = 100;


% 7 dimensional matrix, make some placeholders
x = NaN(4,length(S)); % 6 dimensions as listed in the paper + 1 for LFP voltage

% define initial conditions 
x(1,1) = .5; % half of all receptors bound
x(2,1) = 0.5; % half of all channels open
x(3,1) = 0; % no diffusable factor 
x(4,1) = p.V_rest; % mV??

dt = 1e-3; % hard coded for now

% divide time dependant parameters by the time step
p.ka = p.ka*dt;
p.sa = p.sa*dt;
p.kb = p.kb*dt;
p.a = p.a*dt;
p.b = p.b*dt;

S = S*p.stim_scale;


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

	% now we convert channel opening to voltage using something I made up because this is where the model description ends in the paper
	V = x(4,i-1);
	fx = x(2,i-1)*p.g_OR*(p.V_OR - V) + p.g_leak*(p.V_rest - V);
	x(4,i) =x(4,i-1) + dt*fx;


end	

voltage = x(4,:);
diff_factor = x(3,:);
channels = x(2,:);
receptors = x(1,:);

% now convert this fake LFP into a fake firing rate
t = 1:50; % in ms
K = normpdf(t,10,3).*(-8e-3*t+.0812);
K = K/max(K);
K = K*20;

resp =  (filter(K,1,x(4,:)));

resp = resp + p.f0;


