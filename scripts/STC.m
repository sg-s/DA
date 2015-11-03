% STC.m
% spike-triggered covariance
%
% created by Srinivas Gorur-Shandilya at 3:03 , 30 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [eigen_vectors,eigen_values,C_prior,C_spike] = STC(spikes,S)

% defaults
filter_length = .5; % seconds
dt = 1e-3;

% hash the inputs
temp.filter_length = filter_length;
temp.dt = dt;
temp.spikes = spikes;
temp.S = S;
hash = dataHash(temp);

results = cache(hash);
if ~isempty(results)
	eigen_vectors = results.eigen_vectors;
	eigen_values =  results.eigen_values;
	C_prior = results.C_prior;
	C_spike = results.C_spike;
	return
end

% remove mean from stimulus
S = S-mean(S);

% step 1: compute C_prior
n = floor(filter_length/dt);
C_prior_x = zeros(n,1);
for i = 1:n
	shat = S(n+1:end) - S(n+1 - i:end - i);
	C_prior_x(i) = mean(shat);
end

C_prior = zeros(n);
for i = 1:n
	for j = 1:n
		C_prior(i,j) = C_prior_x(i)*C_prior_x(j);
	end
end

% step 2: now compute C_spike
C_spike_x = zeros(n,1);
spike_times = find(spikes);
spike_times(spike_times<n+1) = [];
for i = 1:n
	shat = S(spike_times) - S(spike_times - i);
	C_spike_x(i) = mean(shat);
end
C_spike_1 = zeros(n);
C_spike_2 = zeros(n);
for i = 1:n
	textbar(i,n)
	shat_x = S(spike_times) - S(spike_times - i);
	for j = 1:n
		shat_y = S(spike_times) - S(spike_times - j);
		C_spike_1(i,j) = mean(shat_x.*shat_y);
		C_spike_2(i,j) = C_spike_x(i)*C_spike_x(j);
	end
end

C_spike = C_spike_1 - C_spike_2;
C = C_spike_1 - C_spike_2 - C_prior;
[eigen_vectors,D] = eigs(C);

eigen_values = (diag(D));
eigen_values = eigen_values/max(eigen_values);

% cache 
results.eigen_vectors = eigen_vectors;
results.eigen_values = eigen_values;
results.C_prior = C_prior;
results.C_spike = C_spike;
cache(hash,[])
cache(hash,results)
