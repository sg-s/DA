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
% p.Kr = @filter_gamma2;
% p.Ka = @filter_gamma2;
% p.Kr_tau1 = 10;
% p.Kr_n = 2;
% p.Kr_tau2 = 20;
% p.Kr_A = 1;
% p.Ka_tau1 = 15;
% p.Ka_n = 2;
% p.Ka_tau2 = 25;
% p.Ka_A = 2;
% p.A = 100; % these are the parameters for the hill function at the output
% p.n = 2;
% p.Kd = 5;
% p.offset = 2;
% % 
% created by Srinivas Gorur-Shandilya at 10:52 , 08 December 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function [f,Sigma,eta] = ReceptorAdaptationModel(s,p)
switch nargin
case 0
	help ReceptorAdaptationModel
	return
case 1
	error('No parameters specified')
end

% revert to default filters if none specified
if ~isfield(p,'Ka')
	p.Ka = @filter_gamma2;
end
if ~isfield(p,'Kr')
	p.Kr = @filter_gamma2;
end

% make the filters
t = 1:1000;
Kr = p.Kr(p.Kr_tau1,p.Kr_n,p.Kr_tau2,p.Kr_A,t);
Ka = p.Ka(p.Ka_tau1,p.Ka_n,p.Ka_tau2,p.Ka_A,t);

% filter the stimulus
s_hat = s - mean(s);

if p.Ka_A < 0
	% code for ignore adaptation filter
	Sigma = s;
else
	Sigma = s./(1 + filter(Ka,1,s_hat));
end

Sigma_hat = Sigma - mean(Sigma);

eta = mean(Sigma) + filter(Kr,1,Sigma_hat);

% pass through output non-linearity
if p.A < 0
	% this is code for ignoring the output nonlinearty
	f = eta;
else
	x = [p.A p.Kd p.n];
	f = hill(x,eta);
end

f = f + p.offset;














