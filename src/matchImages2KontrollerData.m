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

% convert all files to v7.3
disp('Converting files to v7.3...')
for i = 1:length(all_kontroller_files)
	convertMATFileTo73([kontroller_data_path all_kontroller_files(i).name]);
end

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
	image_data.control_roi = cell(1);
	image_data.test_roi = cell(1);
	image_data.image_time = cell(1);
	clc
	load([kontroller_data_path all_kontroller_files(i).name])
	disp([kontroller_data_path all_kontroller_files(i).name])

	% check if the trial info in the timestamps is not fucked
	if length(unique(timestamps(1,:)*100 + timestamps(2,:))) ~= length(timestamps(1,:))
		disp('Timestamps fucked. repairing...')
		unique_paradigms = unique(timestamps(1,:));
		for j = 1:length(unique_paradigms)
			temp= (timestamps(2,timestamps(1,:) ==unique_paradigms(j)));
			if length(temp) > 1
				timestamps(2,timestamps(1,:) ==unique_paradigms(j)) = 1:length(temp);
			end
		end
		disp(timestamps(1:2,:))
	end


	disp('Merging data...')
	for j = 1:size(timestamps,2)
		disp(['Paradigm: ' mat2str(timestamps(1,j))])
		disp(['Trial # : ' mat2str(timestamps(2,j))])
		paradigm = timestamps(1,j);
		trial = timestamps(2,j);


		[time_error,loc] = min(abs(timestamps(3,j) - video_starts));
		time_error = datevec(time_error);
		time_error = time_error(end) + time_error(end-1)*60; % in seconds
		if time_error < 30
			disp('Match found...')

			% load this video file
			clearvars images control_roi test_roi andor_elapsed_time 
			disp('Loading matching video...')
			load([image_path image_files(loc).name]);

			% figure out when the LED (controlled by kontroller) turned on and off
			mean_lux = (squeeze(mean(mean(images,1))));
			d_mean_lux = diff(mean_lux);
			[on_size,led_on] = max(d_mean_lux);
			[off_size,led_off] = min(d_mean_lux);
			if min([on_size -off_size]/std(d_mean_lux)) > 10
				images = images(:,:,led_on:led_off);
				disp('Traces aligned successfully...')

				% calcualte dt for imaging
				try
					image_dt = length(data(paradigm).PID(1,:))*1e-4/size(images,3);
					t = image_dt:image_dt:length(data(paradigm).PID(1,:))*1e-4;

					if exist('control_roi') && exist('test_roi')
						% first make a the masks
						control_roi_mask = sum(control_roi,3);
						test_roi_mask = sum(test_roi,3);
					else
						warning('No Control ROI or test ROI!!')
						control_roi_mask = NaN*images(:,:,1);
						test_roi_mask = NaN*images(:,:,1);
					end

					roi_signals = zeros(size(images,3),1);
					temp = images.*repmat(control_roi_mask,1,1,length(roi_signals));
					temp = (squeeze(sum(sum(temp,1))));

					% divide by the pre-stimulus mean
					temp = temp/mean(temp(3:find(t>5,1,'first')));

					% interpolate to 1ms resolution 
					data(paradigm).GCamp6_control_roi(trial,:) = interp1(t,temp,1e-3:1e-3:max(t));

					% now do the test ROI
					roi_signals = zeros(size(images,3),1);
					temp = images.*repmat(test_roi_mask,1,1,length(roi_signals));
					temp = (squeeze(sum(sum(temp,1))));

					% divide by the pre-stimulus mean
					temp = temp/mean(temp(3:find(t>5,1,'first')));

					% interpolate to 1ms resolution 
					data(paradigm).GCamp6_test_roi(trial,:) = interp1(t,temp,1e-3:1e-3:max(t));
				catch
					warning('Something went horribly wrong, kontroller data is missing!')
				end

			end
		else
			warning('NO IMAGE FOUND THAT MATCHES THIS TRIAL!!!')
			keyboard
		end

	end

	disp('Merging and saving all data...')
	save([kontroller_data_path all_kontroller_files(i).name],'data','-append')
end
