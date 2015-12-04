% matchImages2KontrollerData.m
% finds the image sequence that corresponds to trials in kontroller data, and then matches them and creates a combined mat file. used in calcium imaging experiments, where kontroller runs the experiment and handles ephys data, but a different software (MicroManager) handles the image acquisition 
% 
% created by Srinivas Gorur-Shandilya at 10:43 , 03 December 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


function [] = matchImages2KontrollerData(kontroller_data_path,image_path)

% some input validation
if ~strcmp(kontroller_data_path(end),oss)
	kontroller_data_path = [kontroller_data_path oss];
end
if ~strcmp(image_path(end),oss)
	image_path = [image_path oss];
end


all_kontroller_files = dir([kontroller_data_path,'*.mat']);

% check if all these files are actually kontroller generated
disp('Verifying kontroller-generated data files...')
rm_this = [];
for i = 1:length(all_kontroller_files)
	textbar(i,length(all_kontroller_files))
	variable_names = whos('-file',[kontroller_data_path all_kontroller_files(i).name]);
	if any(strcmp('timestamps',{variable_names.name})) && any(strcmp('data',{variable_names.name})) && any(strcmp('ControlParadigm',{variable_names.name})) 
		% OK
	else
		rm_this = [rm_this i];
	end
end
all_kontroller_files(rm_this) = [];

% first make a matrix of video starts and stops for all the videos we have
disp('Determining when each video starts and stops...')
image_files = dir([image_path, '*.mat']);
video_starts = NaN(length(image_files),1);
video_stops = NaN(length(image_files),1);
for i = 1:length(image_files)
	textbar(i,length(image_files))
	clear absolute_time
	load([image_path image_files(i).name],'absolute_time')
	video_starts(i) = min(absolute_time);
	video_stops(i) = max(absolute_time);
end

% for each of these, load and start matching

for i = 1:length(all_kontroller_files)
	clearvars data ControlParadigm timestamps metadata image_data image_data
	image_data = struct;
	image_data.paradigm = cell(1);
	image_data.trial = cell(1);
	image_data.images = cell(1);
	image_data.control_roi = cell(1);
	image_data.test_roi = cell(1);
	image_data.image_time = cell(1);
	load([kontroller_data_path all_kontroller_files(i).name])
	disp([kontroller_data_path all_kontroller_files(i).name])
	disp('Merging data...')
	for j = 1:size(timestamps,2)
		disp(j)
		[time_error,loc] = min(abs(timestamps(3,j) - video_starts));
		time_error = datevec(time_error);
		time_error = time_error(end) + time_error(end-1)*60; % in seconds
		if time_error < 30

			% load this video file
			clearvars images control_roi test_roi andor_elapsed_time 
			load([image_path image_files(loc).name]);

			% figure out when the LED (controlled by kontroller) turned on and off
			mean_lux = (squeeze(mean(mean(images,1))));
			d_mean_lux = diff(mean_lux);
			[on_size,led_on] = max(d_mean_lux);
			[off_size,led_off] = min(d_mean_lux);
			if min([on_size -off_size]/std(d_mean_lux)) > 10
				image_time = andor_elapsed_time(led_on:led_off);
				image_time = image_time - min(image_time);

				% combine the data
				image_data.paradigm{end+1} = timestamps(1,j);
				image_data.trial{end+1} = timestamps(2,j);
				image_data.images{end+1} = images(:,:,led_on:led_off);
				try
					image_data.control_roi{end+1} = control_roi;
				catch
					image_data.control_roi{end+1} = 0;
					warning('NO CONTROL ROI!!!. IMAGE FILE IS:')
					disp(image_files(loc).name)
				end
				try
					image_data.test_roi{end+1} = test_roi;
				catch
					image_data.test_roi{end+1} = 0;
					warning('NO TEST ROI!!! image file is:')
					disp(image_files(loc).name)
				end
				image_data.image_time{end+1} = image_time;
			end
		else
			warning('NO IMAGE FOUND THAT MATCHES THIS TRIAL!!!')
		end

	end

	% determine if this is v7.3 or not
	convertMATFileTo73([kontroller_data_path all_kontroller_files(i).name]);


	disp('Merging and saving all data...')
	save([kontroller_data_path all_kontroller_files(i).name],'image_data','-append')
end
