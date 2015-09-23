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

% reshape the LFP signals
block_length = 1e4;
reshaped_LFP = LFP(1e4:end-1e4-1,:);
reshaped_LFP = reshape(reshaped_LFP,block_length,width(reshaped_LFP)*length(reshaped_LFP)/block_length);

% also reshape the PID
reshaped_PID = PID(1e4:end-1e4-1,:);
reshaped_PID = reshape(reshaped_PID,block_length,width(reshaped_PID)*length(reshaped_PID)/block_length);

% remove a mean from the first epoch
for i = 1:width(reshaped_LFP)
	reshaped_LFP(:,i) = reshaped_LFP(:,i) - mean(reshaped_LFP(1:5e3,i));
end

% make colour scheme for block analysis
filter_length = 1000;
offset = 200;
all_offsets = [1 3 5.8 7.8];
window_length = 2;

figure('outerposition',[0 0 600 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
c = parula(5);
yy = mean(mean(reshaped_LFP)) + mean(3*std(reshaped_LFP));
for i = 1:length(all_offsets)
	plot([all_offsets(i) all_offsets(i)+window_length],[yy yy],'Color',c(i,:),'LineWidth',10);
end

errorShade(1e-3*(1:length(reshaped_LFP)),mean2(reshaped_LFP),sem(reshaped_LFP));
xlabel('Time since high \rightarrow low switch (s)')
ylabel('\DeltaLFP (mV)')

PrettyFig

if being_published
	snapnow
	delete(gcf)
end


%%
% We now extract filters in two second blocks in this triggered time (starting from the time of switch from high to low). 


% let's try to pull out filters from every epoch
sr = 1e3; % sampling rate, Hz
K = cache(DataHash([reshaped_LFP,reshaped_PID]));
ft = -99:700;
if isempty(K)
	K = NaN(length(all_offsets),filter_length-offset,width(reshaped_LFP));
	for i = 1:width(reshaped_LFP)
		textbar(i,width(reshaped_PID))
		for j = 1:length(all_offsets)

			x = NaN(length(reshaped_PID),1);
			y = x;

			a = 1+ all_offsets(j)*sr - filter_length + offset;
			z = all_offsets(j)*sr + offset + window_length*sr;
			this_stim = reshaped_PID(a:z,i);
			ff = fit((1:length(this_stim))',this_stim,'poly1');
			this_stim = this_stim - ff(1:length(this_stim));
			x(a:z) = this_stim;

			a = all_offsets(j)*sr;
			z = a + window_length*sr;
			this_resp = reshaped_LFP(a:z,i);
			ff = fit((1:length(this_resp))',this_resp,'poly1');
			this_resp = this_resp - ff(1:length(this_resp));
			y(a:z) = this_resp;

			otp = false(length(y),1);
			otp(a:z) = true;

			[this_K,ft] = fitFilter2Data(x,y,'reg',1,'offset',offset,'OnlyThesePoints',otp,'filter_length',filter_length);
			K(j,:,i) = this_K(100:end-102);
			ft = ft(100:end-102);
		end
	end
	cache(DataHash([reshaped_LFP,reshaped_PID]),K);
end


% make linear predictions on the detrended data
LFP_pred = NaN*reshaped_LFP;
offset = 100;
for i = 1:width(reshaped_LFP)
	for j = 1:length(all_offsets)
		a = all_offsets(j)*1e3;
		z = a + window_length*1e3;
		x = reshaped_PID(:,i);
		this_pred = convolve(1e-3*(1:length(x)),x,squeeze(K(j,:,i)),ft);
		LFP_pred(a:z,i) = this_pred(a:z)  - mean(this_pred(a:z));
	end
end

% remove trend from output in order to compare it to the linear prediction. 
for i = 1:width(reshaped_LFP)
	for j = 1:length(all_offsets)
		a = all_offsets(j)*1e3;
		z = a + window_length*1e3;
		ff = fit((1:z-a)',reshaped_LFP(a:z-1,i),'poly1');
		reshaped_LFP(a:z-1,i) = reshaped_LFP(a:z-1,i) - ff(1:z-a);
	end
end

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1), hold on
filtertime = 1e-3*ft;
for i = 1:length(all_offsets)
	errorShade(filtertime,mean2(squeeze(K(i,:,:))),sem(squeeze(K(i,:,:))),'Color',c(i,:));
end
xlabel('Filter Lag (s)')
ylabel('Filter Amplitude')

subplot(1,3,2), hold on
L = {};
r2 = NaN(width(LFP_pred),length(all_offsets));
for j = 1:length(all_offsets)
	a = all_offsets(j)*sr;
	z = window_length*1e3 + a;
	for i = 1:width(LFP_pred)
		r2(i,j) = rsquare(LFP_pred(a:z,i),reshaped_LFP(a:z,i));
	end
	x = j + 0.1*randn(width(reshaped_LFP),1);
	plot(x,r2(:,j),'+','Color',c(j,:));
	L{j} = ['t= ' oval(all_offsets(j)) 's'];
end
set(gca,'XTick',[1:4],'XTickLabel',L)
ylabel('r^2')
title('Linear Prediction Quality')

% plot responses when prediction is good
subplot(1,3,3), hold on
for j = 1:length(all_offsets)
	a = all_offsets(j)*sr;
	z = window_length*1e3 + a;

	x = (LFP_pred(a:z,r2(:,1) > .5));
	y = (reshaped_LFP(a:z,r2(:,1) > .5));
	x = x(:); y = y(:);
	ff = fit(x,y,'poly1');
	plot(linspace(min(x),max(x),30),ff(linspace(min(x),max(x),30)),'Color',c(j,:))
end

xlabel('Linear Prediction')
ylabel('\DeltaLFP (mV)')

PrettyFig()

if being_published
	snapnow
	delete(gcf)
end

% ######## #### ########  #### ##    ##  ######      ########     ###    ######## ########  ######  
% ##        ##  ##     ##  ##  ###   ## ##    ##     ##     ##   ## ##      ##    ##       ##    ## 
% ##        ##  ##     ##  ##  ####  ## ##           ##     ##  ##   ##     ##    ##       ##       
% ######    ##  ########   ##  ## ## ## ##   ####    ########  ##     ##    ##    ######    ######  
% ##        ##  ##   ##    ##  ##  #### ##    ##     ##   ##   #########    ##    ##             ## 
% ##        ##  ##    ##   ##  ##   ### ##    ##     ##    ##  ##     ##    ##    ##       ##    ## 
% ##       #### ##     ## #### ##    ##  ######      ##     ## ##     ##    ##    ########  ######  


%% Firing Rate Analysis
% We now perform a similar analysis, but for the firing rates. 

% reshape the firing rate signals
block_length = 1e4;
reshaped_fA = fA(1e4:end-1e4-1,:);
reshaped_fA = reshape(reshaped_fA,block_length,width(reshaped_fA)*length(reshaped_fA)/block_length);


% remove a mean from the first epoch
for i = 1:width(reshaped_fA)
	reshaped_fA(:,i) = reshaped_fA(:,i) - mean(reshaped_fA(1:5e3,i));
end

% make colour scheme for block analysis
filter_length = 1000;
offset = 200;
all_offsets = [1 3 5.8 7.8];
window_length = 2;

figure('outerposition',[0 0 600 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
c = parula(5);
yy = mean(mean(reshaped_fA)) + mean(3*std(reshaped_fA));
for i = 1:length(all_offsets)
	plot([all_offsets(i) all_offsets(i)+window_length],[yy yy],'Color',c(i,:),'LineWidth',10);
end

errorShade(1e-3*(1:length(reshaped_fA)),mean2(reshaped_fA),sem(reshaped_fA));
xlabel('Time since high \rightarrow low switch (s)')
ylabel('Firing Rate (\DeltaHz)')

PrettyFig

if being_published
	snapnow
	delete(gcf)
end


%%
% We now extract filters in two second blocks in this triggered time (starting from the time of switch from high to low). 


% let's try to pull out filters from every epoch
sr = 1e3; % sampling rate, Hz
K2 = cache(DataHash([reshaped_fA,reshaped_PID]));
if isempty(K2)
	K2 = NaN(length(all_offsets),filter_length-offset,width(reshaped_fA));
	for i = 1:width(reshaped_fA)
		textbar(i,width(reshaped_PID))
		for j = 1:length(all_offsets)

			x = NaN(length(reshaped_PID),1);
			y = x;

			a = 1+ all_offsets(j)*sr - filter_length + offset;
			z = all_offsets(j)*sr + offset + window_length*sr;
			this_stim = reshaped_PID(a:z,i);
			ff = fit((1:length(this_stim))',this_stim,'poly1');
			this_stim = this_stim - ff(1:length(this_stim));
			x(a:z) = this_stim;

			a = all_offsets(j)*sr;
			z = a + window_length*sr;
			this_resp = reshaped_fA(a:z,i);
			ff = fit((1:length(this_resp))',this_resp,'poly1');
			this_resp = this_resp - ff(1:length(this_resp));
			y(a:z) = this_resp;

			otp = false(length(y),1);
			otp(a:z) = true;

			[this_K2,ft] = fitFilter2Data(x,y,'reg',1,'offset',offset,'OnlyThesePoints',otp,'filter_length',filter_length);
			K2(j,:,i) = this_K2(100:end-102);
			ft = ft(100:end-102);
		end
	end
	cache(DataHash([reshaped_fA,reshaped_PID]),K2);
end


% make linear predictions on the detrended data
fp = NaN*reshaped_fA;
offset = 100;
for i = 1:width(reshaped_fA)
	for j = 1:length(all_offsets)
		a = all_offsets(j)*1e3;
		z = a + window_length*1e3;
		x = reshaped_PID(:,i);
		this_pred = convolve(1e-3*(1:length(x)),x,squeeze(K2(j,:,i)),ft);
		fp(a:z,i) = this_pred(a:z)  - mean(this_pred(a:z));
	end
end

% remove trend from output in order to compare it to the linear prediction. 
for i = 1:width(reshaped_fA)
	for j = 1:length(all_offsets)
		a = all_offsets(j)*1e3;
		z = a + window_length*1e3;
		ff = fit((1:z-a)',reshaped_fA(a:z-1,i),'poly1');
		reshaped_fA(a:z-1,i) = reshaped_fA(a:z-1,i) - ff(1:z-a);
	end
end

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1), hold on
filtertime = 1e-3*ft;
for i = 1:length(all_offsets)
	errorShade(filtertime,mean2(squeeze(K2(i,:,:))),sem(squeeze(K2(i,:,:))),'Color',c(i,:));
end
xlabel('Filter Lag (s)')
ylabel('Filter Amplitude')

subplot(1,3,2), hold on
L = {};
r2 = NaN(width(fp),length(all_offsets));
for j = 1:length(all_offsets)
	a = all_offsets(j)*sr;
	z = window_length*1e3 + a;
	for i = 1:width(fp)
		r2(i,j) = rsquare(fp(a:z,i),reshaped_fA(a:z,i));
	end
	x = j + 0.1*randn(width(reshaped_fA),1);
	plot(x,r2(:,j),'+','Color',c(j,:));
	L{j} = ['t= ' oval(all_offsets(j)) 's'];
end
set(gca,'XTick',[1:4],'XTickLabel',L)
ylabel('r^2')
title('Linear Prediction Quality')

% plot responses when prediction is good
subplot(1,3,3), hold on
for j = 1:length(all_offsets)
	a = all_offsets(j)*sr;
	z = window_length*1e3 + a;

	x = (fp(a:z,r2(:,1) > .5));
	y = (reshaped_fA(a:z,r2(:,1) > .5));
	x = x(:); y = y(:);
	ff = fit(x,y,'poly1');
	plot(linspace(min(x),max(x),30),ff(linspace(min(x),max(x),30)),'Color',c(j,:))
end

xlabel('Linear Prediction')
ylabel('Firing Rate (\DeltaHz)')

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
