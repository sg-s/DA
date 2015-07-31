% condSpikeProb.m
% spiking probability, conditional on spike
% note that this is not the same as the ISI distribution
% 
% created by Srinivas Gorur-Shandilya at 11:05 , 31 July 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [p_spike,bin_centres] = condSpikeProb(A,T_max,dt)

if ~isvector(A)
	error('needs a vector to work with')
end

T_max = floor(T_max/1e-4);
dt = floor(dt/1e-4);

A = A(:);
spiketimes = find(A);

% make the x axis
bin_left = 0:dt:T_max;
bin_right =bin_left + dt;
y = zeros(length(bin_left),1);

for i = 1:length(spiketimes)
	spike_offsets = spiketimes - spiketimes(i);
	spike_offsets(spike_offsets<1) = [];
	spike_offsets(spike_offsets>T_max) = [];
	for j = 1:length(spike_offsets)
		y(find(spike_offsets(j)>bin_left,1,'last')) = y(find(spike_offsets(j)>bin_left,1,'last')) + 1;
	end

end
p_spike = y/length(spiketimes);
bin_centres = (bin_left + bin_right)/2;
bin_centres = bin_centres*1e-4;

% lose the last point
p_spike(end) = [];
bin_centres(end) = [];