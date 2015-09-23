% LiuWangReduced.m
% a reduced version of the Liu-Wang model with 6 parameters:
% 
% p.tau_m 		-- membrane time constant
% p.g 			-- dimensionless number, ratio of g_AHP/g_L
% p.r_k 		-- r_k = V_K - V_rest (voltage)
% p.tau_Ca 		-- calcium buffering time constant
% p.r_reset  	-- reset voltage (after emitting a spike)
% p.Alpha 		-- calcium influx / spike
% p.r_th 		-- threshold for spiking
% 
% in addition, there are some additional parameters governing how the stimulus is converted into a current:
% 
% p.A  			-- a simple scaling on the stimulus
%
% created by Srinivas Gorur-Shandilya at 10:06 , 23 April 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function [f,spikes,r,Ca] = LiuWangReduced(stimulus,p)
if ~nargin
	help LiuWangReduced
	return
end

% constrains
lb.g = 0;
lb.Alpha = 0;
ub.r_th = 0;
ub.r_k = 0;
ub.r_reset = 0;
lb.tau_m = 0;
lb.tau_Ca = 0;

% relative constraints
if p.tau_Ca < p.tau_m
	temp = p.tau_m;
	p.tau_m = p.tau_Ca;
	p.tau_Ca = temp;
end


time = 1e-3*(1:length(stimulus));
stimulus = stimulus(:);
spikes = 0*time;
r = 0*time;
Ca = 1*time;

% initial conditions
r(1) = p.r_reset;
dt = (time(2) - time(1));

for i = 2:length(time)
	% convert stimulus into a current
	I_stim = stimulus(i-1)*p.A;

	% calculate derivative
	d = -r(i-1) + I_stim - p.g*(Ca(i-1))*(r(i-1) - p.r_k);

	% update variables using the Euler approximation
	r(i) = r(i-1) + dt*d/p.tau_m;
	Ca(i) = Ca(i-1)*(1- dt/p.tau_Ca);

	if Ca(i) < 0
		Ca(i)=0; % Ca can't be negative!
	end

	% handle spiking
	if r(i) >= p.r_th
		if r(i-1) < 0 && ~spikes(i-1)
			r(i-1) = 30; % fake a spike. this is purely cosmetic
			r(i) = p.r_reset;
			spikes(i)=1;
			Ca(i) = Ca(i) + p.Alpha; % amount of calcium entering/spike
		else
			r(i) = p.r_reset;
		end
	end
end

% generate the firing rate estimate 
f = spiketimes2f(spikes,time,1e-3);
