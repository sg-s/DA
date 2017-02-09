% contrastSteepnessModel
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

function R = contrastSteepnessModel(S,p)

% list parameters for clarity
p.k0;
p.B;
p.n;
p.tau;

% bounds
lb.k0 = 19; % tweak this based on what the maximum steepness can be
lb.B = 0;
lb.n = 1;
lb.tau = 1;

ub.tau = 1e3;


filter_length = 4*(p.n*p.tau);
t = 0:filter_length; 
Kg = generate_simple_filter(p.tau,p.n,t);
K = NaN*S;

%% variant 1 ----  K x [S]+   ------------------------------------------------------
S(S<0) = 0;
for i = 1:width(S)
	Shat = filter(Kg,1,S(:,i));
	K(:,i) = p.k0./(1 + p.B*Shat);
end

%% variant 2 ----  -K x [S]-  ------------------------------------------------------
% S(S>0) = 0; S = -S;
% for i = 1:width(S)
% 	Shat = filter(Kg,1,S(:,i));
% 	K(:,i) = p.k0./(1 + p.B*Shat);
% end

% %% variant 3 ----  [K x S]+   ------------------------------------------------------
% for i = 1:width(S)
% 	Shat = filter(Kg,1,S(:,i));
% 	Shat(Shat<0) = 0;
% 	K(:,i) = p.k0./(1 + p.B*Shat);
% end

% %% variant 4 ----  -[K x S]-  ------------------------------------------------------
% for i = 1:width(S)
% 	Shat = filter(Kg,1,S(:,i));
% 	Shat(Shat>0) = 0; Shat = -Shat;
% 	K(:,i) = p.k0./(1 + p.B*Shat);
% end


R(1) = nanmean(nanmean(K(1e3:4e3,:)));
R(2) = nanmean(nanmean(K(6e3:9e3,:)));
R(3) = R(1) - R(2);

function f = generate_simple_filter(tau,n,t)
f = t.^n.*exp(-t/tau); % functional form in paper
f = f/tau^(n+1)/gamma(n+1); % normalize appropriately






 	