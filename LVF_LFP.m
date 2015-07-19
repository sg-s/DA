% LVF_LFP.m
% recreates fig3 with LFP recordings
% 
% created by Srinivas Gorur-Shandilya at 11:52 , 14 July 2015. Contact me at http://srinivas.gs/contact/
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

%% 
% In this document, we attempt to recreate the results from Fig 3 showing fast gain changes with a widely distributed stimulus. 

%% Stimulus 
% First, we check that the stimulus is the same as the one we used six months ago (in LargeVarianceFlickering.m). The following plot shows the old stimulus (blue) and the new stimulus (red). 

% load old stimulus
load('/local-data/DA-paper/large-variance-flicker/2015_01_28_CS_ab3_2_EA.mat')
PID = data(4).PID;

load('/local-data/DA-paper/large-variance-flicker/2015_01_28_CS_ab3_3_EA.mat')
PID = vertcat(PID,data(4).PID);
PID = PID';
PID = PID(1:10:end,:);
% remove baseline
baseline = mean(PID(1:300,1));
PID = PID - baseline;

% load new stimulus
load('/local-data/DA-paper/large-variance-flicker/2015_07_14_LVFtest.mat')
newPID = data(4).PID';
newPID = newPID(1:10:end,:);
time = 1e-3*(1:length(PID));

% remove baseline
baseline = mean(newPID(1:300,1));
newPID = newPID - baseline;

% plot
clear l
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
l(1) = errorShade(time(1:10:end),mean2(PID(1:10:end,:)),sem(PID(1:10:end,:)),'Color',[0 0 1]);
l(2) = errorShade(time(1:10:end),mean2(newPID(1:10:end,:)),sem(newPID(1:10:end,:)),'Color',[1 0 0]);
xlabel('Time (s)')
ylabel('PID (V)')
set(gca,'XLim',[20 60])
legend(l,{'January 28 2015','July 14 2015'})

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%%
% It looks very similar, but smaller, which is expected, as the PID gets less sensitive over time. What if we rescale it? 

ff = fit(mean2(newPID),mean2(PID),'poly1');

clear l
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
l(1) = errorShade(time(1:10:end),mean2(PID(1:10:end,:)),sem(PID(1:10:end,:)),'Color',[0 0 1]);
l(2) = errorShade(time(1:10:end),ff(mean2(newPID(1:10:end,:))),sem(newPID(1:10:end,:)),'Color',[1 0 0]);
xlabel('Time (s)')
ylabel('PID (V)')
legend(l,{'January 28 2015','July 14 2015 Rescaled'})
set(gca,'XLim',[20 60])
PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%% Example Data
% How does the neuron respond to this stimulus? Does the LFP and/or firing rates show evidence for fast gain control? In the following figure, we plot the stimulus, LFP and the response from all the neurons in the dataset. The shading shows the standard error of the mean. Each neuron is in a different colour. 

p = '/local-data/DA-paper/large-variance-flicker/LFP/';
[PID, LFP, fA, paradigm, orn, AllControlParadigms, paradigm_hashes,sequence] = consolidateData(p,1);

% set to NaN firing rates that are 0
fA(:,max(fA) == 0) = NaN;

% throw out trials where we didn't record the LFP, for whatever reason
not_LFP = 0*orn;
for i = 1:width(LFP)
	not_LFP(i) = abs(mean2(LFP(:,i)));
end
LFP(:,not_LFP< 0.5) = NaN;

% filter the LFP
filteredLFP = LFP;
for i = 1:width(LFP)
	a = find(~isnan(LFP(:,i)),1,'first');
	z = find(~isnan(LFP(:,i)),1,'last');
	if isempty(a)
		a = 1;
	end
	if isempty(z)
		z = length(LFP);
	end
	filteredLFP(a:z,i) = 10*filter_trace(LFP(a:z,i),1000,10);
end

figure('outerposition',[0 0 900 700],'PaperUnits','points','PaperSize',[900 700]); hold on
clear a
for i = 1:3
	a(i) = subplot(3,1,i); hold on
	set(a(i),'XLim',[20 50])
end
time = 1e-3*(1:length(PID));
c = parula(length(unique(orn))+1);
for i = 1:length(unique(orn))
	example_orn = i;
	axes(a(1))
	errorShade(time,mean2(PID(:,orn==example_orn)),sem(PID(:,orn==example_orn)),'Color',c(i,:),'SubSample',50);
	axes(a(2))
	errorShade(time,mean2(filteredLFP(:,orn==example_orn)),sem(filteredLFP(:,orn==example_orn)),'Color',c(i,:),'SubSample',50);
	axes(a(3))
	errorShade(time,mean2(fA(:,orn==example_orn)),sem(fA(:,orn==example_orn)),'Color',c(i,:),'SubSample',50);
	
end
ylabel(a(1),'Stimulus (V)')
ylabel(a(2),'\DeltaLFP (mV)')
ylabel(a(3),'Firing Rate (Hz)')
 set(a(2),'YLim',[-4 2.5])

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%% Example Filters
% We now extract all possible filters from this dataset, for each neuron: 

% K -- PID -> LFP filter
if ~exist('K','var')
	K = NaN(1e3,length(unique(orn)));
	for i = 1:length(unique(orn))
		resp = mean2(filteredLFP(20e3:55e3,orn==i));
		stim = mean2(PID(20e3:55e3,orn==i));
		stim(1:400) = [];
		resp(end-399:end) = [];
		temp = FitFilter2Data(stim,resp,[],'reg=1;','filter_length=1399;');
		% throw out 200ms on either end
		temp(1:200) = [];
		temp(end-199:end) = [];
		K(:,i) = temp;
	end
end

if ~exist('K2','var')
	K2 = NaN(1e3,length(unique(orn)));
	for i = 1:length(unique(orn))
		stim = mean2(filteredLFP(20e3:55e3,orn==i));
		resp = mean2(fA(20e3:55e3,orn==i));
		stim(1:400) = [];
		resp(end-399:end) = [];
		temp = FitFilter2Data(stim,resp,[],'reg=1;','filter_length=1399;');
		% throw out 200ms on either end
		temp(1:200) = [];
		temp(end-199:end) = [];
		K2(:,i) = temp;
	end
end

if ~exist('K3','var')
	K3 = NaN(1e3,length(unique(orn)));
	for i = 1:length(unique(orn))
		stim = mean2(PID(20e3:55e3,orn==i));
		resp = mean2(fA(20e3:55e3,orn==i));
		stim(1:400) = [];
		resp(end-399:end) = [];
		temp = FitFilter2Data(stim,resp,[],'reg=1;','filter_length=1399;');
		% throw out 200ms on either end
		temp(1:200) = [];
		temp(end-199:end) = [];
		K3(:,i) = temp;
	end
end

figure('outerposition',[0 0 1400 450],'PaperUnits','points','PaperSize',[1400 550]); hold on
clear a
for i = 1:3
	a(i) = subplot(1,3,i); hold on
	xlabel('Lag (s)')
end
filtertime = 1e-3*(1:length(K)) - .2;
for i = 1:width(K)
	axes(a(1))
	plot(filtertime,K)
	axes(a(2))
	plot(filtertime,K2)
	axes(a(3))
	plot(filtertime,K3)

end
title(a(1),'PID \rightarrow LFP Filter')
title(a(2),'LFP \rightarrow Firing rate Filter')
title(a(3),'PID \rightarrow Firing rate Filter')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%% Example Data: Fast Gain Control
% First, we check if we see evidence for fast gain control in the firing rate output of the neuron. 



% make LN predictions of the LFP and the response
fp = NaN(length(fA),length(unique(orn)));
for i = 1:length(unique(orn))
	filtertime = 1e-3*(1:1e3)-.2;
	fp(:,i) = convolve(time,mean2(PID(:,orn==i)),K3(:,i),filtertime);
	% correct for some trivial scaling
	a = fp(:,i);
	b = mean2(fA(:,orn == i));
	temp  = isnan(a) | isnan(b);
	temp = fit(a(~temp),b(~temp),'poly1'); 
	fp(:,i) = temp(fp(:,i));
end

LFP_pred = NaN(length(fA),length(unique(orn)));
for i = 1:length(unique(orn))
	filtertime = 1e-3*(1:1e3)-.2;
	LFP_pred(:,i) = convolve(time,mean2(PID(:,orn==i)),K(:,i),filtertime);
	% correct for some trivial scaling
	a = LFP_pred(:,i);
	b = mean2(filteredLFP(:,orn == i));
	temp  = isnan(a) | isnan(b);
	temp = fit(a(~temp),b(~temp),'poly7'); 
	LFP_pred(:,i) = temp(LFP_pred(:,i));
end

for i = 1:length(unique(orn))

	% gain analysis -- LN model
	figure('outerposition',[0 0 700 700],'PaperUnits','points','PaperSize',[1000 700]); hold on

	% first we do the firing rates
	clear ph
	ph(3) = subplot(2,2,1); hold on
	ph(4) = subplot(2,2,2); hold on

	hl_min = .1;
	hl_max = 10;
	history_lengths = logspace(log10(hl_min),log10(hl_max),30);

	resp = mean2(fA(:,orn==i));
	stim = mean2(PID(:,orn==i));
	pred = (fp(:,i));

	stim = stim(20e3:55e3);
	pred = pred(20e3:55e3);
	resp = resp(20e3:55e3);

	history_lengths = findValidHistoryLengths(1e-3,stim,pred,resp,30,.33);

	[p,~,~,~,~,history_lengths]=GainAnalysisWrapper2('response',resp,'prediction',pred,'stimulus',stim,'time',1e-3*(1:length(resp)),'ph',ph,'history_lengths',history_lengths,'use_cache',1,'engine',@GainAnalysis5,'example_history_length',history_lengths(11));
	title(ph(4),['ORN ',mat2str(i),' Firing Rate'])

	% now do the same for the LFP
	ph(3) = subplot(2,2,3); hold on
	ph(4) = subplot(2,2,4); hold on

	resp = mean2(filteredLFP(:,orn==i));
	pred = LFP_pred(:,i);
	pred = pred(20e3:55e3);
	resp = resp(20e3:55e3);

	[p,~,~,~,~,history_lengths]=GainAnalysisWrapper2('response',resp,'prediction',pred,'stimulus',stim,'time',1e-3*(1:length(resp)),'ph',ph,'history_lengths',history_lengths,'example_history_length',history_lengths(11),'use_cache',1,'engine',@GainAnalysis5);
	title(ph(4),['ORN ',mat2str(i),' LFP'])
	xlabel(ph(3),'LFP (pred)')
	ylabel(ph(3),'LFP (data)')

	PrettyFig;

	if being_published
		snapnow
		delete(gcf)
	end
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
