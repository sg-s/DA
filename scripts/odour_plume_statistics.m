
pHeader;

%% Odour plume statistics
% In this document, I look at the statistics of real odour plumes. I generate them using a fan blowing over a vial of either Apple Cider Vinegar or ethyl acetate. 

% get data 
root = '/data/DA-paper/data-for-paper/odor-plumes/';
root = '~/Desktop/odor-plumes/';

% first, gather the data
allfiles = dir([root '*.kontroller']);
allfiles(cellfun(@(x) strcmp(x(1),'.'), {allfiles.name})) = [];

odour_x_location = [];
odour_y_location = [];
fan_power = [];
PID_location = [];
odour_name = {};
PID = [];


for i = 1:length(allfiles)
	clear data spikes ControlParadigm

	load([root allfiles(i).name],'-mat')
	%disp({ControlParadigm.Name})

	a = strfind(allfiles(i).name,'PID_')+4;
	z = strfind(allfiles(i).name,'cm')-1;
	this_PID_location = str2double(allfiles(i).name(a:z));

	this_PID = vertcat(data.PID);
	this_PID = this_PID';
	this_PID = this_PID(1:10:end,:);

	z = strfind(allfiles(i).name,'_fan')-1;
	this_odour = allfiles(i).name(12:z);

	this_odour = repmat({this_odour},size(this_PID,2),1);
	this_PID_location = repmat(this_PID_location,size(this_PID,2),1);

	% figure out the positions of the odour for each trial
	n_trials_per_paradigm = cellfun(@width,{data.PID});
	this_x = [];
	this_y = [];
	this_fan_power = [];
	for j = 1:length(data)
		disp(ControlParadigm(j).Name)
		z = strfind(ControlParadigm(j).Name,'cm'); 
		if length(z) > 1
			yz = z(2);
			z = z(1);
			a = strfind(ControlParadigm(j).Name,'_'); ya = a(end); a = a(1);
			x = str2double(ControlParadigm(j).Name(a+1:z-1));
			y = str2double(ControlParadigm(j).Name(ya+1:yz-1));
			this_x = [this_x; zeros(n_trials_per_paradigm(j),1)+x];
			this_y = [this_y; zeros(n_trials_per_paradigm(j),1)+y];
		else
			a = strfind(ControlParadigm(j).Name,'_');
			x = str2double(ControlParadigm(j).Name(a+1:z-1));
			this_x = [this_x; zeros(n_trials_per_paradigm(j),1)+x];
			this_y = [this_y; zeros(n_trials_per_paradigm(j),1)];
		end

		% figure out fan power
		if isempty(strfind(ControlParadigm(j).Name,'fan_off'))
			this_fan_power = [this_fan_power; ones(n_trials_per_paradigm(j),1)];
		else
			this_fan_power = [this_fan_power; zeros(n_trials_per_paradigm(j),1)];
		end
		
	end


	% consolidate
	PID = [PID this_PID];
	PID_location = [PID_location; this_PID_location];
	odour_x_location = [odour_x_location; this_x];
	odour_y_location = [odour_y_location; this_y];
	fan_power = [fan_power; this_fan_power];
	odour_name = [odour_name; this_odour];

end

clearvars -except PID odour_name PID_location odour_x_location odour_y_location being_published fan_power

nplots = unique(sum([odour_x_location 1e3*odour_y_location],2));


prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end


%% Version Info
%
pFooter;


