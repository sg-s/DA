% analyseCalciumSignals
% works on the output of matchImages2KontrollerData, and converts the complicated calcium images into simple intensity measurements
% 
% created by Srinivas Gorur-Shandilya at 5:14 , 03 December 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [] = analyseCalciumSignals(path_name)

% some input validation
if ~strcmp(path_name(end),oss)
	path_name = [path_name oss];
end

allfiles = dir([path_name,'*.mat']);

for i = 1:length(allfiles)
	disp('Loading...')
	disp([path_name allfiles(i).name])
	load([path_name allfiles(i).name])
	disp('DONE. Computing the calcium fluorescence...')

	for paradigm = 1:length(data)
		for trial = 1:width(data(paradigm).PID)
			disp(['Paradigm: ' oval(paradigm) 'Trial:' oval(trial)])

			% find the corresponding entry in image_data
			use_this = [];
			for j = 1:length(image_data.paradigm)
				if image_data.paradigm{j} == paradigm & image_data.trial{j} == trial
					use_this = j;
				end
			end
			
			if isempty(use_this)
				warning('No image data!!!')

			else
				% now for each control ROI, build a vector 
				control_roi_mask = image_data.control_roi{use_this};
				roi_signals = zeros(size(image_data.images{use_this},3),size(control_roi_mask,3));

				% calcualte dt for imaging
				image_dt = length(data(paradigm).PID(trial,:))*1e-4/length(roi_signals);

				% assume odour turns on @ 5s
				t = image_dt:image_dt:length(data(paradigm).PID(trial,:))*1e-4;

				for j = 1:width(roi_signals)
					temp = image_data.images{use_this}.*repmat(control_roi_mask(:,:,j),1,1,length(roi_signals));
					roi_signals(:,j) = (squeeze(sum(sum(temp,1))));

					% divide by the pre-stimulus mean
					roi_signals(:,j) = roi_signals(:,j)/mean(roi_signals(3:find(t>5,1,'first'),j));
				end
				data(paradigm).GCamp_control(trial,:) = mean(roi_signals,2);

				test_roi_mask = image_data.test_roi{use_this};
				roi_signals = zeros(size(image_data.images{use_this},3),size(test_roi_mask,3));
				for j = 1:width(roi_signals)
					temp = image_data.images{use_this}.*repmat(test_roi_mask(:,:,j),1,1,length(roi_signals));
					roi_signals(:,j) = (squeeze(sum(sum(temp,1))));

					% divide by the pre-stimulus mean
					roi_signals(:,j) = roi_signals(:,j)/mean(roi_signals(3:find(t>5,1,'first'),j));
				end

				data(paradigm).GCamp(trial,:) = mean(roi_signals,2);
				data(paradigm).GCamp_time(trial,:) = t;

			end
		end
	end

	save([path_name allfiles(i).name],'data','-append')

end
