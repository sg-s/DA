% contrastSteepnessModel2
% just like contrastSteepnessModel, but has a stimulus scale and offset term
% 
% this model accepts a bunch of input signals, that carry information about the contrast of the signal
% typically, this would be the derivative of the LFP (or something like this)
% it also assumes that the signal is 10s long, with the switch occuring at 5s
% it uses this predictor signal to calculate K, the dynamic steepness. it then averages this, 
% and returns a two-element vector (k1,k2)
% which is the steepness in the two epochs (high and low variance)
% 
% created by Srinivas Gorur-Shandilya at 1:49 , 16 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function R = contrastSteepnessModel2(S,p)

% list parameters for clarity
p.k0;
p.B;
p.n;
p.tau;

p.s0;
p.s1;

% bounds
lb.k0 = 0;
lb.B = 0;
lb.n = 1;
lb.tau = 1;

ub.tau = 500;
ub.n = 2;


filter_length = 4*(p.n*p.tau);
t = 0:filter_length; 
Kg = generate_simple_filter(p.tau,p.n,t);

S = S*p.s0 + p.s1;

K = NaN*S;
S(S<0) = 0;

for i = 1:width(S)
	% do hi stimulus
	Shat = filter(Kg,1,S(:,i));
	K(:,i) = p.k0./(1 + p.B*Shat);
end

R(1) = nanmean(nanmean(K(1e3:5e3,:)));
R(2) = nanmean(nanmean(K(6e3:9e3,:)));

function f = generate_simple_filter(tau,n,t)
f = t.^n.*exp(-t/tau); % functional form in paper
f = f/tau^(n+1)/gamma(n+1); % normalize appropriately






 	