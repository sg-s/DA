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

	% remove a quadratic trend from the stimulus
	stimulus = data(i).PID;
	time = 1:length(stimulus);
	for j = 1:width(stimulus)
		temp = fit(time(:),stimulus(:,j),'poly2');
		stimulus(:,j) = stimulus(:,j) - temp(time) + mean(stimulus(:,j));
	end

	orn_data(i).stimulus = stimulus;
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
do_these = [7 9 13 15 16 7 8 10 14 17 18 11 17 25];
for i = 1:length(do_these)
	do_this = do_these(i);
	disp(do_this)
	orn_data(do_this) = computeInstGain(orn_data(do_this));
end