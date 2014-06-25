% XJWNeuronEuler.m
% part of the DA project
% simulates and solves using a first order Euler method Liu and Wang's neuron model.
% parameter structure:
% p.A -- scaling constant for stimulus --> current transformation 
% p.Vrest -- resting potential of neuron
% p.gL -- leak conductance
% p.Vreset -- reset potential
% p.gAHP -- AHP (K+) conductance 
% p.Vk -- K channel potential 
% p.tau_Ca -- time constant for calcium buffering 
% p.C -- calcium kick for each spike
% p.Cm -- membrance capcaitance 
% p.Vth -- threshold voltage
% 
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function [spikes] = XJWNeuronEuler(time,stimulus,p)


time = time(:);
stimulus = stimulus(:);
spikes = 0*time;
V = 0*time;
Ca = 1*time;

% initial conditions
V(1) = p.Vrest;


for i = 2:length(time)

	% calculate derivative
	f1 = -p.gL*(V(i-1) - p.Vreset);
	f2 = p.Cm*p.A*stimulus(i-1);
	f3 = -p.gAHP*(Ca(i-1))*(V(i-1) - p.Vk);

	f = f1+f2+f3;

	dt = (time(i) - time(i-1));

	V(i) = V(i-1) + dt*f/p.Cm;
	Ca(i) = Ca(i-1)*(1- dt/p.tau_Ca);

	if V(i) > p.Vth
		V(i-1) = 30; % fake a spike
		V(i) = p.Vreset;
		spikes(i)=1;
		Ca(i) = Ca(i) + p.C;
	end
end


