% temporary script to convert carlotta's data into the new ORNData class
% 
% created by Srinivas Gorur-Shandilya at 1:54 , 26 January 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

% load and convert

combined_data_file = ('/local-data/DA-paper/carlotta-martelli/flickering-stim/data.mat');
load(combined_data_file)

orn_data = ORNData;

for i = 1:length(data)
	textbar(i,length(data))
	orn_data(i).firing_rate = data(i).fA;
	orn_data(i).stimulus = data(i).PID;
	orn_data(i).valve = data(i).Valve;
	orn_data(i).spikes = data(i).spikes;
	orn_data(i).neuron_name = data(i).neuron_name;
	orn_data(i).odour_name = data(i).odour_name;
	orn_data(i).original_name = data(i).original_name;
	orn_data(i).data_creator = 'Carlotta Martelli';
end

% compute filters
for i = 1:length(orn_data)
	orn_data(i) = backOutFilters(orn_data(i));
end

% compute inst. gain
for i = 1:length(orn_data)
	orn_data(i) = computeInstGain(orn_data(i));
end