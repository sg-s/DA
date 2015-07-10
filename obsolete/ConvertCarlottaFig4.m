% ConvertCarlottaFig4
% 
% created by Srinivas Gorur-Shandilya at 3:59 , 10 March 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
%
% this script converts carlotta's fig4 data for 2ac into a more usable format

allfiles= dir('/local-data/DA-paper/carlotta-martelli/fig4/*SPK*.mat');

clear data
for i = 1:length(allfiles)
	% load data
	load(strcat('/local-data/DA-paper/carlotta-martelli/fig4/',allfiles(i).name))
	load(strcat('/local-data/DA-paper/carlotta-martelli/fig4/',strrep(allfiles(i).name,'_SPK','')))

	disp(allfiles(i).name)

	data(i).name = allfiles(i).name;
	data(i).dilutions = dilutions;

	fA = [];
	time = 1e-4*(1:length(PID));
	fA = NaN(length(dilutions),size(tSPKA,2),length(time)/10);
	t = 1e-3*(1:length(fA));
	thisPID = fA;
	for j = 1:length(dilutions)
		textbar(j,length(dilutions))
		for k = 1:size(tSPKA,2)
			fA(j,k,:) = spiketimes2f(squeeze(tSPKA(j,k,:)),time);
			thisPID(j,k,:) = interp1(time,squeeze(PID(j,k,:)),t);
		end
	end

	% subtract the baseline from the PID
	for j = 1:size(tSPKA,2)
		thisPID(:,j,:) = thisPID(:,j,:) -  mean(PID_back(j,1:4e4));
	end

	data(i).fA = fA;
	data(i).PID = thisPID;

end

% now compute metrics for each trial
a = 1050;
z = 1500; % nominal stimulus start and stop
background_stim = NaN(1e4,1);
foreground_stim = NaN(1e4,1);
ORN_peak = NaN(1e4,1);
ORN_peak_time = NaN(1e4,1);
stim_half_time = NaN(1e4,1); % time it takes to go to half max
resp_half_time = NaN(1e4,1); % time it takes to go to half max
ORN_half_time = NaN(1e4,1);
c = 1;


% convert data into PulseData for other uses
PulseData = struct;
for i = 1:length(dilutions)
	PulseData(i).resp = ((squeeze(data(1).fA(i,:,:)))');
	PulseData(i).stim = ((squeeze(data(1).PID(i,:,:)))');
end

% check for bad data
for i = 1:length(PulseData)
	rm_this_trial = [];
	for j = 1:width(PulseData(i).resp)
		if std(PulseData(i).resp(:,j))< 1e-1
			rm_this_trial = [rm_this_trial j];
		end
	end
	PulseData(i).resp(:,rm_this_trial) = [];
	PulseData(i).stim(:,rm_this_trial) = [];
end

for i = 1:length(PulseData)
	rm_this_trial = [];
	for j = 1:width(PulseData(i).resp)
		if std(PulseData(i).stim(:,j)) == 0
			rm_this_trial = [rm_this_trial j];
		end
	end
	PulseData(i).resp(:,rm_this_trial) = [];
	PulseData(i).stim(:,rm_this_trial) = [];
end

save('PulseData.mat','PulseData')



for i = 1:length(data)
	for j = 1:length(data(i).dilutions)
		for k = 1:size(data(i).fA,2)
			if std(squeeze(data(i).PID(j,k,:))) > 1e-6
				if i > 1
					background_stim(c) = mean(data(i).PID(j,k,1:a-100));
				else
					background_stim(c) = 0;
				end
				foreground_stim(c) = mean(data(i).PID(j,k,a+200:z));

				[ORN_peak(c),ORN_peak_time(c)]= max(squeeze((data(i).fA(j,k,a:z))));
				if foreground_stim(c) > background_stim(c)
					this_stim = squeeze(data(i).PID(j,k,:));
					this_stim = this_stim - mean(this_stim(1:a));
					this_stim = this_stim/max(this_stim(a:z));
					this_resp = squeeze(data(i).fA(j,k,:));
					this_resp = this_resp -  mean(this_resp(1:a));
					this_resp = this_resp/max(this_resp(a:z));
					stim_half_time(c) = max([find(this_stim(a:z)>.5,1,'first') NaN]);
					resp_half_time(c) = max([find(this_resp(a:z)>.5,1,'first') NaN]);
				end
				c=c+1;
			end
		end
	end
end
background_stim(c:end) = [];
foreground_stim(c:end) = [];
ORN_peak(c:end) = [];
ORN_peak_time(c:end) = [];
stim_half_time(c:end) = [];
resp_half_time(c:end) = [];

save('2ac_timing.mat','background_stim','foreground_stim','ORN_peak','ORN_peak_time','stim_half_time','resp_half_time')
