% AssembleCMFlickeringData
% 
% this script assembles all the flickering stimulus data from carlotta's paper into a more usable data structure. 
% this data structure has the following sub-fields:
% data.original_name
% data.PID -- stimulus, matrix (trial x time)
% data.fA -- response, matrix
% data.neuron_name
% data.odour_name
% data.K -- may be a matrix, in which case it is the trial-wise filters
% data.LinearFit 
% data.dt -- usually 1e-3 s
% data.history_lengths -- which history lengths should we use for this data set?
%
% created by Srinivas Gorur-Shandilya at 11:28 , 12 March 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
%


datasource = '/local-data/DA-paper/carlotta-martelli/flickering-stim/raw/*.mat';
allfiles = dir(datasource);
consolidated_data_location = '/local-data/DA-paper/carlotta-martelli/flickering-stim/data.mat';

neuron_names = {'ab3A','pb1A','ab3a','pb1a'};
odor_names = {'i5ac','2but','5ol','1but','2ac','d2succ','1o3ol'};

if exist(consolidated_data_location)
	load(consolidated_data_location)
else
	data = struct;
	data.original_name = [];
	data.PID = [];%-- stimulus, matrix (trial x time)
	data.fA = []; %-- response, matrix
	data.neuron_name = [];
	data.odour_name = [];
	data.K = []; %-- may be a matrix, in which case it is the trial-wise filters
	data.LinearFit = [];
	data.dt = []; 
	data.history_lengths = [];
	data.Valve = [];
	data.remove = []; % this stores the remove array in the original data.
end

for i = 1:length(allfiles)
	% fine load it
	filename=(strcat(datasource(1:end-5),allfiles(i).name));
	disp('Prepping:')
	disp(filename)
	clear remove PID stim_signal stA
	load(filename)
	% remove(remove>width(PID)) = [];
	PID = squeeze(PID);
	% PID(remove,:) =[];
	stA = squeeze(stA);
	baddata = (sum(stA') == 0);

	% throw out bad data -- no spikes.
	stA(baddata,:) = [];
	PID(baddata,:) = [];

	% stA(remove,:) = [];
	stim_signal = squeeze(stim_signal);
	if ~isvector(stim_signal)
		stim_signal = mean2(stim_signal);
	end
	stim_signal(stim_signal>.5) = 1;
	stim_signal(stim_signal<1)= 0;
	time = 1e-4*(1:length(PID));

	% clip data
	a = find(stim_signal>0,1,'first');
	a = round(a/1e4)*1e3;
	z = find(stim_signal>0,1,'last');
	z = round(z/1e4)*1e3;

	data(i).fA = spiketimes2f(squeeze(stA),time);
	data(i).fA = data(i).fA(a:z,:);
	data(i).original_name = allfiles(i).name;
	data(i).dt = 1e-3;

	% subsample PID
	if size(PID,2) > size(PID,1)
		PID  = PID';
	end
	data(i).PID = NaN*data(i).fA;
	t = a*1e-3:data(i).dt:z*1e-3;
	for j = 1:width(PID)
		data(i).PID(:,j) = interp1(time,PID(:,j),t);
		data(i).PID(:,j) = data(i).PID(:,j) - mean(mean(PID(1:a*10,j)));
	end


	data(i).Valve = stim_signal(a:z);


	% figure out the odor
	this_odor = [];
	for k = 1:length(odor_names)
		if any(strfind(allfiles(i).name,odor_names{k}))
			this_odor = odor_names{k};
		end
	end
	clear k

	% figure out the neuron
	this_neuron = [];
	for k = 1:length(neuron_names)
		if any(strfind(allfiles(i).name,neuron_names{k}))
			this_neuron = neuron_names{k};
		end
	end
	clear k

	data(i).neuron_name=this_neuron;
	data(i).odour_name = this_odor;
	data(i).remove = remove;
end


% clean up some data
data(11).fA(:,end) = [];
data(11).PID(:,end) =[];
data(17).fA(:,end-1) = [];
data(17).PID(:,end-1) = [];

save(consolidated_data_location,'data')
