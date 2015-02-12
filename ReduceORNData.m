% ReduceORNData.m
% scans a folder for Kontroller files, and automatically combines them all, creating data structures as needed.
% assumes all the data files are driven by the same control paradigms (identical in all aspects)
% created by Srinivas Gorur-Shandilya at 8:52 , 13 November 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function combined_data = ReduceORNData(data_root,allfiles)

dt = 1e-3; % the dt of the output data that this function will return

% % get all the variables to be combined
% for i = 1:length(allfiles)
% 	load(strcat(data_root,allfiles(i).name))
% 	variable_names = fieldnames(data);
% 	for j = 1:length(variable_names)
% 		if ~exist(variable_names{i})
% 			eval(strcat(variable_names{i},'=[];'))
% 		end
% 	end
% end	

fA = [];
fB = [];
PID = [];
neuron = [];
paradigm = {};

% now combine all the data
for i = 1:length(allfiles)
	disp(allfiles(i).name)
	clear data
	clear spikes
	load(strcat(data_root,allfiles(i).name))
	if exist('spikes')
		if any(find(strcmp('A', fieldnames(spikes))))
			% ok, has spike data
			for j = 1:length(spikes)
				for k = 1:width(spikes(j).A)
					disp([j k])
					if length(spikes(j).A(k,:)) > 2
						if sum(spikes(j).A(k,:)) > 1
							% haz spikes
							time = 1e-4*(1:length(data(j).PID(k,:)));
							n = length(data(j).PID(k,:));
							fA = [fA  spiketimes2f(spikes(j).A(k,1:n),time)];
							fB = [fB  spiketimes2f(spikes(j).B(k,1:n),time)];

							paradigm = [paradigm ControlParadigm(j).Name];

							PID = [PID ; interp1(time,data(j).PID(k,:),dt*(1:length(fA)))];

							neuron = [neuron; i];



						end
					end
				end

			end
		end
	end
end

combined_data.fA = fA;
combined_data.fB = fB;
combined_data.PID = PID;
combined_data.paradigm = paradigm;
combined_data.neuron = neuron;
