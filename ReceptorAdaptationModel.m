% ReceptorAdaptationModel.m
% this is a model of ORN response where the ORN response f is given by
% f = N(Rx(s/(s+Kxs)))
% where
% s is the stimulus
% K is a filter that governs how the receptors adapt
% R is a filter that governs the impulse response function of the neuron
% N is a hill function
% and x is meant to indicate convolution
% K and R can be any parametrised filter
% 
% usage:
% f = ReceptorAdaptationModel(s,p)
% where s is a stimulus vector
% and p is a structure with the parameters of this model defined e.g.:
% p.tau1= 8.4062;
% p.K_n= 2.2969;
% p.tau2= 54.3438;
% p.K_A= 1;
% p.a_tau1= 15;
% p.   a_n= 2;
% p.a_tau2= 25;
% p.   a_A= 1;
% p.beta = -1;
% p.      A= 1;
% p.      n= 2;
% p.     Kd= 1;
% p. offset= 1;
% % 
% created by Srinivas Gorur-Shandilya at 10:52 , 08 December 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function [f,shat,shat2] = ReceptorAdaptationModel(s,p)
switch nargin
case 0
	help ReceptorAdaptationModel
	return
case 1
	error('No parameters specified')
end

% make the filters
t = 1:300;
Kr = filter_gamma2(p.tau1,p.K_n,p.tau2,p.K_A,t);
Ka = filter_gamma2(p.a_tau1,p.a_n,p.a_tau2,p.a_A,t);

% transform input into the adapted input
shat = s./(1+p.beta*filter(Ka,1,s-mean(s)));


% run through response filter
shat2 = filter(Kr,1,shat-mean(shat));

% add a offset
shat2 = shat2 + p.offset;

% pass through the non-linearity
if p.n > 0
	x = [p.A p.Kd p.n];
	f = hill(x,shat2);
else
	% only scaling, no nonlinearity
	f = shat2*p.A;
end








