% NWModelEuler.m
% Euler-integration implementation of Nagel-Wilson Model
% 
% created by Srinivas Gorur-Shandilya at 6:24 , 15 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [R,x] = NWModelEuler(S,p)

% bounds
lb.ka = eps;
lb.theta = eps;
lb.a = eps;
lb.b = eps;
lb.sa = eps;
lb.sb = eps;


% 7 dimensional matrix, make some placeholders
x = NaN(7,length(S)); % 6 dimensions as listed in the paper + 1 for LFP voltage

% define initial conditions 
x(1,1) = 1;
x(2:4,1) = 0; % this sets the inactive receptor population to maximum, everything else to 0
x(5,1) = 0; % all channels closed to start with
x(6,1) = 1; % no diffusable factor 
x(7,1) = -10; % mV??

dt = 1e-3; % hard coded for now

% divide time dependant parameters by the time step
p.ka = p.ka/dt;
p.sa = p.sa/dt;
p.kb = p.kb/dt;
p.sb = p.sb/dt;
p.a = p.a/dt;
p.b = p.b/dt;

S  = S*p.stim_scale;


for i = 2:length(S)
	O = S(i-1); % odor stimulus

	% get the receptor activities right
	fx = [-(p.ka*p.sa+O*p.kb*p.sb)	p.sa 						p.sb  						0;
	p.ka*p.sa 					-(p.sa+O*p.theta*p.kb*p.sb)	0 							p.sb;
	O*p.kb*p.sb 				0     						-(p.sb*p.theta*p.ka*p.sa) 	p.sa;
	0 							O*p.theta*p.kb*p.sb			p.theta*p.ka*p.sa 			-(p.sa+p.sb)];

	% calculate receptor variables
	x(1:4,i) = x(1:4,i-1) + dt*(fx*x(1:4,i-1)); % simple Euler integration

	% force their sum to be 1, and for all of them to be positive
	x(x(1:4,i)<0,i) = 0;

	if sum(x(1:4,i)) > 1
		x(1:4,i) = x(1:4,i)./(sum(x(1:4,i)));
	end

	% now the channel opening 
	fx = (x(2,i-1)+x(4,i-1))*p.ko*(1/x(6,i-1))*(1-x(5,i-1)) - p.kc*x(5,i-1); % according to Nagel and Wilson
	x(5,i) = x(5,i-1) + dt*fx;

	% now the diffusable factor
	fx = p.a*x(5,i-1) - p.b*x(6,i-1);  
	x(6,i) = x(6,i-1) + dt*fx;

	% now we convert channel opening to voltage using something I made up because this is where the model description ends in the paper
	fx = (p.V0 + p.R*x(5,i-1)*(p.Ec-x(7,i-1)) - x(7,i-1))/p.tau;
	x(7,i) =x(7,i-1) + dt*fx;

end	

% now convert this fake LFP into a fake firing rate
t = 1:50; % in ms
K = normpdf(t,10,3).*(-8e-3*t+.0812);
K = K/max(K);
K = K*20;

R =  (filter(K,1,x(7,:)));
R = R + p.f0;

