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
% How does the neuron respond to this stimulus? Does the LFP and/or firing rates show evidence for fast gain control? In the following figure, we plot the stimulus, LFP and the response from one example neuron. The shading shows the standard error of the mean:

p = '/local-data/DA-paper/large-variance-flicker/LFP/';
[PID, LFP, fA, paradigm, orn, AllControlParadigms, paradigm_hashes,sequence] = consolidateData(p,1);

% sort the paradigms sensibly
sort_value = [];
for i = 1:length(AllControlParadigms)
	sort_value(i) = mean(mean(AllControlParadigms(i).Outputs));
end

[~,idx] = sort(sort_value);
AllControlParadigms = AllControlParadigms(idx);
paradigm_new = paradigm*NaN;
for i = 1:length(idx)
	paradigm_new(paradigm == idx(i)) = i;
end
paradigm = paradigm_new;

% throw out trials where we didn't record the LFP, for whatever reason
not_LFP = 0*orn;
for i = 1:width(LFP)
	not_LFP(i) = abs(mean2(LFP(:,i)));
end
LFP(:,not_LFP< 0.5) = NaN;

% filter the LFP
filteredLFP = LFP;
for i = 1:width(LFP)
	filteredLFP(:,i) = 10*filter_trace(LFP(:,i),1000,10);
end

figure('outerposition',[0 0 900 700],'PaperUnits','points','PaperSize',[900 700]); hold on
subplot(3,1,1), hold on
time = 1e-3*(1:length(PID));
errorShade(time,mean2(PID(:,orn==1)),sem(PID(:,orn==1)),'Color',[0 0 0]);
set(gca,'XLim',[20 50])
ylabel('Stimulus (V)')

subplot(3,1,2), hold on
errorShade(time,mean2(filteredLFP(:,orn==1)),sem(filteredLFP(:,orn==1)),'Color',[0 0 0]);
set(gca,'XLim',[20 50])
ylabel('LFP (mV)')

subplot(3,1,3), hold on
errorShade(time,mean2(fA(:,orn==1)),sem(fA(:,orn==1)),'Color',[0 0 0]);
set(gca,'XLim',[20 50])
ylabel('Firing Rate (Hz)')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%% Example Filters
% We now extract all possible filters from this dataset:

if ~exist('K','var')
	K = NaN(1e3,length(orn));
	for i = 1:length(orn)
		resp = filteredLFP(20e3:55e3,i);
		rm_this = isnan(resp);
		resp(rm_this) = [];
		if length(resp) 
			stim = PID(20e3:55e3,i);
			stim(rm_this) = [];
			stim(1:400) = [];
			resp(end-399:end) = [];
			K(:,i) = FitFilter2Data(stim,resp,[],'reg=1;','filter_length=999;');
		end
	end
end

if ~exist('K2','var')
	K2 = NaN(1e3,length(orn));
	for i = 1:length(orn)
		resp = fA(20e3:55e3,i);
		stim = filteredLFP(20e3:55e3,i);
		rm_this = isnan(resp) | isnan(stim);
		resp(rm_this) = [];
		stim(rm_this) = [];
		if length(resp) 
			stim = filter_trace(stim,1e3,10);
			stim(1:400) = [];
			resp(end-399:end) = [];
			K2(:,i) = FitFilter2Data(stim,resp,[],'reg=1;','filter_length=999;');
		end
	end
end

if ~exist('K3','var')
	K3 = NaN(1e3,length(orn));
	for i = 1:length(orn)
		resp = fA(20e3:55e3,i);
		stim = PID(20e3:55e3,i);
		rm_this = isnan(resp) | isnan(stim);
		resp(rm_this) = [];
		stim(rm_this) = [];
		if length(resp) & sum(resp) > 0
			stim(1:400) = [];
			resp(end-399:end) = [];
			K3(:,i) = FitFilter2Data(stim,resp,[],'reg=1;','filter_length=999;');
		end
	end
end

figure('outerposition',[0 0 1400 450],'PaperUnits','points','PaperSize',[1400 450]); hold on

subplot(1,3,1), hold on
title('PID \rightarrow LFP Filter')
time = 1e-3*(1:501)-.1;
plot_this = find(orn == 1);
plot_this = setdiff(plot_this,find(isnan(sum(K))));
errorShade(time,mean2(K(300:800,plot_this)),sem(K(300:800,plot_this)));
xlabel('Lag (s)')

subplot(1,3,2), hold on
title('LFP \rightarrow Firing rate Filter')
plot_this = find(orn == 1);
plot_this = setdiff(plot_this,find(isnan(sum(K2))));
errorShade(time,mean2(K2(300:800,plot_this)),sem(K2(300:800,plot_this)));
xlabel('Lag (s)')

subplot(1,3,3), hold on
title('PID \rightarrow Firing rate Filter')
plot_this = find(orn == 1);
plot_this = setdiff(plot_this,find(isnan(sum(K3))));
errorShade(time,mean2(K3(300:800,plot_this)),sem(K3(300:800,plot_this)));
xlabel('Lag (s)')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%% Example Data: Fast Gain Control
% First, we check if we see evidence for fast gain control in the firing rate output of the neuron. 


% make linear predictions of the LFP and the response
fp = fA*NaN;
for i = 1:length(orn)
	filtertime = 1e-3*(1:501)-.1;
	this_K = K3(300:800,i);
	fp(:,i) = convolve(time,PID(:,i),this_K,filtertime);
	% correct for some trivial scaling
	a = fp(:,i);
	b = fA(:,i);
	temp  = isnan(a) | isnan(b);
	temp = fit(a(~temp),b(~temp),'poly1'); 
	fp(:,i) = temp(fp(:,i));
end

LFP_pred = filteredLFP*NaN;
for i = 1:length(orn)
	filtertime = 1e-3*(1:501)-.1;
	this_K = K(300:800,i);
	if ~any(isnan(this_K))
		LFP_pred(:,i) = convolve(time,PID(:,i),this_K,filtertime);
		% correct for some trivial scaling
		a = LFP_pred(:,i);
		b = filteredLFP(:,i);
		temp  = isnan(a) | isnan(b);
		temp = fit(a(~temp),b(~temp),'poly1'); 
		LFP_pred(:,i) = temp(LFP_pred(:,i));
	end
end

% gain analysis -- LN model
figure('outerposition',[0 0 700 700],'PaperUnits','points','PaperSize',[1000 700]); hold on

% first we do the firing rates
clear ph
ph(3) = subplot(2,2,1); hold on
ph(4) = subplot(2,2,2); hold on

hl_min = .1;
hl_max = 10;
history_lengths = logspace(log10(hl_min),log10(hl_max),30);

resp = mean2(fA(:,orn==1));
stim = mean2(PID(:,orn==1,:));
pred = mean2(fp(:,orn==1));
stim = stim(20e3:55e3);
pred = pred(20e3:55e3);
resp = resp(20e3:55e3);

temp = fit(pred(:),resp(:),'poly5');
pred = temp(pred);

[p,~,~,~,~,history_lengths]=GainAnalysisWrapper2('response',resp,'prediction',pred,'stimulus',stim,'time',1e-3*(1:length(resp)),'ph',ph,'history_lengths',history_lengths,'example_history_length',history_lengths(7),'use_cache',1,'engine',@GainAnalysis5);
title(ph(4),'ORN 1 Firing Rate')

% now do the same for the LFP
ph(3) = subplot(2,2,3); hold on
ph(4) = subplot(2,2,4); hold on

resp = mean2(filteredLFP(:,orn==1));
pred = mean2(LFP_pred(:,orn==1));
pred = pred(20e3:55e3);
resp = resp(20e3:55e3);

temp = fit(pred(:),resp(:),'poly5');
pred = temp(pred);

[p,~,~,~,~,history_lengths]=GainAnalysisWrapper2('response',resp,'prediction',pred,'stimulus',stim,'time',1e-3*(1:length(resp)),'ph',ph,'history_lengths',history_lengths,'example_history_length',history_lengths(7),'use_cache',1,'engine',@GainAnalysis5);
title(ph(4),'ORN 1 LFP')

PrettyFig;

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
