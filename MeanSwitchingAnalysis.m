% MeanSwitchingAnalysis.m
% analysis of mean switching experiment
% 
% created by Srinivas Gorur-Shandilya at 2:34 , 27 August 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


% this code determines if this function is being called by publish() or not
calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end
tic

%% Mean Switching Analysis
% In this document, we analyse the responses of ORNs to stimuli where we rapidly switch from two ensembles of stimuli, differing only in their means. 

%% Stimulus
% In the following figure, we show what the stimulus looks like: 

p = '/local-data/switching/mean/use-this/';
[PID, LFP, fA, paradigm, orn] = consolidateData(p,1);

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,4,1:3), hold on
time = 1e-3*(1:length(PID));
plot(time,PID(:,1),'k');
set(gca,'XLim',[30 100])
xlabel('Time (s)')
ylabel('Stimulus (V)')

subplot(1,4,4), hold on
x = 0:0.01:2.5;
y = NaN(length(x),width(PID));
a = 10e3+1:10e3:length(PID)-5e3;
for i = 1:length(a)
	this_start = a(i);
	try
		temp = PID(this_start-5e3:this_start,2);
		y(:,i) = hist(temp(1:10:end),x);
		y(:,i) = y(:,i)/sum(y(:,i));
	end
end
errorShade(x,mean2(y),sem(y),'Color',[1 0 0],'Shading',0.4);

x = 0:0.01:2.5;
y = NaN(length(x),width(PID));
a = 15e3+1:10e3:length(PID)-5e3;
for i = 1:length(a)
	this_start = a(i);
	try
		temp = PID(this_start-4e3:this_start,2);
		y(:,i) = hist(temp(1:10:end),x);
		y(:,i) = y(:,i)/sum(y(:,i));
	end
end
errorShade(x,mean2(y),sem(y),'Color',[0 0 1],'Shading',0.4);
xlabel('Stimulus (V)')
ylabel('Probability')

PrettyFig()

if being_published
	snapnow
	delete(gcf)
end


%% LFP Analysis
% We now line up all the LFP signals by the switching time, and look at how the LFP changes in aggregate when we switch from a low mean stimulus to a high mean stimulus. Since the LFP drifts over time, and we don't really care about that, we construct a triggered LFP where we subtract the mean of the LFP during the low stimulus window in each epoch and then average. 

% remove spikes from LFP
for i = 1:width(LFP)
	LFP(:,i) = filtfilt(ones(20,1),20,LFP(:,i))*10; % in mV
end



triggeredLFP = NaN(1e4,0);
a = 20e3;
z = length(LFP)-5e3;
for start_here = a:10e3:z
	for i = 1:width(LFP)
		this_segment = LFP(start_here:start_here+1e4-1,i);
		this_segment = this_segment - mean(this_segment(1:5e3));
		triggeredLFP = [triggeredLFP this_segment];
	end
end

% make colour scheme for block analysis
all_offsets = [1 3 6 8];
window_length = 2;

figure('outerposition',[0 0 600 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
c = parula(5);
yy = mean(mean(triggeredLFP)) + mean(2*std(triggeredLFP));
for i = 1:length(all_offsets)
	plot([all_offsets(i) all_offsets(i)+window_length],[yy yy],'Color',c(i,:),'LineWidth',10);
end

errorShade(1e-3*(1:length(triggeredLFP)),mean2(triggeredLFP),sem(triggeredLFP));
xlabel('Time since high \rightarrow low switch (s)')
ylabel('\DeltaLFP (mV)')

PrettyFig()

if being_published
	snapnow
	delete(gcf)
end

%%
% We now extract filters in two second blocks in this triggered time (starting from the time of switch from high to low). 


K = NaN(length(all_offsets),1e3,width(PID));
filter_length = 1200;
offset = 200;
for i = 1:width(PID)
	for j = 1:length(all_offsets)
		x = PID(:,i);
		y = LFP(:,i);
		OnlyThesePoints = zeros(length(PID),1);
		a = 20e3 + all_offsets(j)*1e3;
		z = length(LFP)-10e3;
		for start_here = a:10e3:z
			OnlyThesePoints(start_here:start_here+window_length*1e3) = 1;

			this_stim = x(start_here-filter_length+offset:start_here-filter_length+offset+2e3);
			this_stim = this_stim - mean(this_stim);
			this_stim = this_stim/std(this_stim);
			% ff = fit((1:length(this_stim))',this_stim,'poly1');
			% this_stim = this_stim - ff(1:length(this_stim));
			x(start_here-filter_length+offset:start_here-filter_length+offset+window_length*1e3) = this_stim;

			this_resp = y(start_here:start_here+window_length*1e3);
			this_resp = this_resp - mean(this_resp);
			this_resp = this_resp/std(this_resp);
			% ff = fit((1:length(this_resp))',this_resp,'poly1');
			% this_resp = this_resp - ff(1:length(this_resp));
			y(start_here:start_here+window_length*1e3) = this_resp;
		end

		y(~OnlyThesePoints) =  NaN;

		[this_K,ft] = fitFilter2Data(x,y,'reg',1,'offset',offset,'OnlyThesePoints',logical(OnlyThesePoints),'filter_length',filter_length,'normalise',false);
		K(j,:,i) = this_K(100:end-102);
		ft = ft(100:end-102);
	end
end

figure('outerposition',[0 0 600 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
filtertime = 1e-3*ft;
for i = 1:length(all_offsets)
	errorShade(filtertime,mean2(squeeze(K(i,:,:))),sem(squeeze(K(i,:,:))),'Color',c(i,:));
end
xlabel('Filter Lag (s)')
ylabel('Filter Amplitude')

PrettyFig()

if being_published
	snapnow
	delete(gcf)
end

%% Firing Rate Analysis
% We now perform a similar analysis, but for the firing rates. 


triggeredfA = NaN(1e4,0);
a = 20e3;
z = length(LFP)-5e3;
for start_here = a:10e3:z
	for i = 1:width(LFP)
		this_segment = fA(start_here:start_here+1e4-1,i);
		this_segment = this_segment - mean(this_segment(1:5e3));
		triggeredfA = [triggeredfA this_segment];
	end
end

% make colour scheme for block analysis
all_offsets = [1 3 6 8];
window_length = 2;

figure('outerposition',[0 0 600 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
c = parula(5);
yy = mean(mean(triggeredfA)) + mean(2*std(triggeredfA));
for i = 1:length(all_offsets)
	plot([all_offsets(i) all_offsets(i)+window_length],[yy yy],'Color',c(i,:),'LineWidth',10);
end

errorShade(1e-3*(1:length(triggeredfA)),mean2(triggeredfA),sem(triggeredfA));
xlabel('Time since high \rightarrow low switch (s)')
ylabel('\Delta Firing Rate (Hz)')

PrettyFig()

if being_published
	snapnow
	delete(gcf)
end

%%
% We now extract LN models in two second blocks in this triggered time (starting from the time of switch from high to low). On the left are filters, and on the right are linear fits to the residuals. 


K2 = NaN(length(all_offsets),1e3,width(PID));
filter_length = 1200;
offset = 200;
for i = 1:width(PID)
	for j = 1:length(all_offsets)
		x = PID(:,i);
		y = fA(:,i);
		OnlyThesePoints = zeros(length(PID),1);
		a = 20e3 + all_offsets(j)*1e3;
		z = length(fA)-10e3;
		for start_here = a:10e3:z
			OnlyThesePoints(start_here:start_here+window_length*1e3) = 1;

			this_stim = x(start_here-filter_length+offset:start_here-filter_length+offset+2e3);
			this_stim = this_stim - mean(this_stim);
			this_stim = this_stim/std(this_stim);
			% ff = fit((1:length(this_stim))',this_stim,'poly1');
			% this_stim = this_stim - ff(1:length(this_stim));
			x(start_here-filter_length+offset:start_here-filter_length+offset+window_length*1e3) = this_stim;

			this_resp = y(start_here:start_here+window_length*1e3);
			this_resp = this_resp - mean(this_resp);
			this_resp = this_resp/std(this_resp);
			% ff = fit((1:length(this_resp))',this_resp,'poly1');
			% this_resp = this_resp - ff(1:length(this_resp));
			y(start_here:start_here+window_length*1e3) = this_resp;
		end

		y(~OnlyThesePoints) =  NaN;

		[this_K,ft] = fitFilter2Data(x,y,'reg',1,'offset',offset,'OnlyThesePoints',logical(OnlyThesePoints),'filter_length',filter_length,'normalise',false);
		K2(j,:,i) = this_K(100:end-102);
		ft = ft(100:end-102);
	end
end

% make linear predictions 
fp = NaN*fA;
offset = 100;
for i = 1:width(fA)
	for j = 1:length(all_offsets)
		x = PID(:,i);
		y = fA(:,i);
		y = y - mean(y);
		y = y/std(y);
		a = 20e3 + all_offsets(j)*1e3;
		z = length(fA)-10e3;
		this_pred = convolve(1e-3*(1:length(x)),x,squeeze(K2(j,:,i)),ft);
		for start_here = a:10e3:z
			fp(start_here:start_here+window_length*1e3,i) = this_pred(start_here:start_here+window_length*1e3);
		end
	end
end

% plot fp vs fA only for the first block
reshaped_fp = reshape(fp(:,1),1e4,length(fp)/1e4);
reshaped_fA = reshape(fA(:,1),1e4,length(fp)/1e4);
reshaped_fp(:,1:2) = []; reshaped_fA(:,1:2) = [];
reshaped_fp(:,end) = []; reshaped_fA(:,end) = [];


figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
filtertime = 1e-3*ft;
for i = 1:length(all_offsets)
	errorShade(filtertime,mean2(squeeze(K2(i,:,:))),sem(squeeze(K2(i,:,:))),'Color',c(i,:));
end
xlabel('Filter Lag (s)')
ylabel('Filter Amplitude')

subplot(1,2,2), hold on
for i = 1:length(all_offsets)
	temp_fp = reshaped_fp(all_offsets(i)*1e3:(all_offsets(i)+window_length)*1e3,:);
	temp_fp = temp_fp(:);
	temp_fA = reshaped_fA(all_offsets(i)*1e3:(all_offsets(i)+window_length)*1e3,:);
	temp_fA = temp_fA(:);
	ff = fit(temp_fp,temp_fA,'poly1');
	plot([min(temp_fp) max(temp_fp)],ff([min(temp_fp) max(temp_fp)]),'Color',c(i,:))
end
xlabel('Linear Prediction')
ylabel('Neuron Firing Rate (Hz)')


PrettyFig()

if being_published
	snapnow
	delete(gcf)
end



%% Version Info
% The file that generated this document is called:
disp(mfilename)

%%
% and its md5 hash is:
Opt.Input = 'file';
disp(DataHash(strcat(mfilename,'.m'),Opt))

%%
% This file should be in this commit:
[status,m]=unix('git rev-parse HEAD');
if ~status
	disp(m)
end

t = toc;

%% 
% This document was built in: 
disp(strcat(oval(t,3),' seconds.'))

% tag the file as being published 
% add homebrew path
path1 = getenv('PATH');
path1 = [path1 ':/usr/local/bin'];
setenv('PATH', path1);

if being_published
	unix(strjoin({'tag -a published',which(mfilename)}));
end
