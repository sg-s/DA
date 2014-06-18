% AssembleData
% assembles data from different experiments into one massive structure that DAPaper will then use
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

datasource{1} = '/data/orn/data-for-da-paper';
datasource{2} = '/data/orn/data-for-da-paper/flickering-stim-carlotta';

neuron_names = {'ab3A','pb1A','ab3a','pb1a'};
odor_names = {'i5ac','2but','5ol','1but','2ac','d2succ','1o3ol'};

consolidated_data_location = '/local-data/DA-paper/data.mat';
if exist('/local-data/DA-paper/data.mat')
	load('/local-data/DA-paper/data.mat')
end

for i = 1:length(datasource)
	allfiles = dir(datasource{i});
	for j = 1:length(allfiles)
		if isempty(strmatch(allfiles(j).name(1),'.')) && ~allfiles(j).isdir
			if isempty(cell2mat(strfind({data.original_name},allfiles(j).name)))
				% fine load it
				filename=(strcat(datasource{i},'/',allfiles(j).name));
				disp('Prepping:')
				disp(filename)
				[PID, time, f, Valve, uncropped] = PrepData3(filename);
				data(end+1).original_name = allfiles(j).name;

				% figure out the odor
				this_odor = [];
				for k = 1:length(odor_names)
					if any(strfind('final_2011_05_16_ab3A_5ol3X-3_20ml_30sec_30ms_rand.mat',odor_names{k}))
						this_odor = odor_names{k};
					end
				end
				clear k

				% figure out the neuron
				this_neuron = [];
				for k = 1:length(neuron_names)
					if any(strfind('final_2011_05_16_ab3A_5ol3X-3_20ml_30sec_30ms_rand.mat',neuron_names{k}))
						this_neuron = neuron_names{k};
					end
				end
				clear k

				data(end).neuron=this_neuron;
				data(end).odor = this_odor;
				data(end).PID = PID;
				data(end).time = time;
				data(end).ORN = f;
				data(end).Valve=  Valve;
				data(end).full_data = uncropped;


			end
		end
	end

end
clear i

save('/local-data/DA-paper/data.mat','data','-append')