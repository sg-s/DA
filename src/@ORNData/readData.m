% ORNData/readData.m
% readData only works when all the paradigms have the same length
% this is really meant to work with just one paradigm
% 
% created by Srinivas Gorur-Shandilya at 9:03 , 16 February 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [od] = readData(od,path_name,use_cache)

assert(exist(path_name,'dir')==7,'First argument should be a valid path')

% check the path
if strcmp(path_name(end),oss)
else
	path_name = [path_name oss];
end

allfiles = dir([path_name '*.mat']);

% remove some junk
rm_this = [];
for i = 1:length(allfiles)
	variable_names = whos('-file',[path_name allfiles(i).name]);
	if any(strcmp('spikes',{variable_names.name})) && any(strcmp('data',{variable_names.name})) && any(strcmp('ControlParadigm',{variable_names.name}))
	else
		rm_this = [rm_this i];
	end
end
allfiles(rm_this) = [];

for i = 1:length(allfiles)
	od(i) = ORNData;
	disp([path_name allfiles(i).name])
	load([path_name allfiles(i).name])

	for j = 1:length(data)
		if ~isempty(data(j).PID)
			for k = 1:width(data(j).PID)
				% load the stimulus
				this_pid = data(j).PID(k,1:10:end); 
				od(i).stimulus(:,width(od(i).stimulus)+1) = this_pid(:);
				% load the LFP
				this_LFP = data(j).voltage(k,1:10:end); this_LFP = this_LFP(:);
				try
					this_LFP = this_LFP - fastFiltFilt(ones(1e4,1),1e4,this_LFP);
					this_LFP = this_LFP*10; % to get the units right, now in mV
				catch
					disp('error in filtering')
					keyboard
				end
				od(i).LFP(:,width(od(i).LFP)+1) = this_LFP;

				% check if we can load spikes
				if width(spikes(j).A) < k
					od(i).firing_rate(:,width(od(i).firing_rate)+1) = this_LFP*NaN;
				else
				% OK
					this_spikes = spikes(j).A(k,:);
					od(i).firing_rate(:,width(od(i).firing_rate)+1) = spiketimes2f(this_spikes,1e-4*(1:length(this_spikes)),1e-3,3e-2);
				end
			end
		end
	end

	% automatically annotate the data so we only use good trials
	od(i).use_these_trials= ~((any(isnan(od(i).LFP))) | (any(isnan(od(i).firing_rate))) | (any(isnan(od(i).stimulus))));

end