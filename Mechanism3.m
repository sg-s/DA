% Mechanism3.m
% in this document, we attempt to understand the mechanism behind gain adaptation
% this is a continuation of experiments detailed in Mechanism2.m 
% 
% created by Srinivas Gorur-Shandilya at 1:15 , 06 April 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end
tic

%% Higher ATR doses works
% Titrating ATR doses and feeding durations leads to a sweetspot where the ORNs respond normally to odour and can be easily activated by light. It is very hit or miss for now, and the control paradigms for each neuron had to be hand-tuned to get it right. 

%% Responses to Odour Flicker with Light Backgrounds
% In the following figures, we present flickering odour stimuli to and record their responses to the odour with and without an additional light activation. Each figure corresponds to data from one neuron. 

datapath = '/local-data/DA-paper/reachr/2015_04_16_RR_F1_ab3_1_EA_1mM.mat';
haz_data = [16:19];

PID = [];
fA = [];
ParadigmNames = {};
paradigm = [];

[stim,resp,ParadigmNames,paradigm] = MechanismAnalysis_PrepData(datapath,haz_data,[],[],{},[],0);
alldata(1).stim = stim;
alldata(1).resp = resp;
alldata(1).ParadigmNames = ParadigmNames;
alldata(1).paradigm = paradigm;
MechanismAnalysis_PlotGain(stim,resp,ParadigmNames,paradigm,0);


PrettyFig();



if being_published
	snapnow
	delete(gcf)
end


%%
% Here's another neuron:

datapath = '/local-data/DA-paper/reachr/2015_04_16_RR_F2_ab3_1_EA_1mM.mat';
haz_data = [16 18];

PID = [];
fA = [];
ParadigmNames = {};
paradigm = [];

[stim,resp,ParadigmNames,paradigm] = MechanismAnalysis_PrepData(datapath,haz_data,[],[],{},[],0);
datapath =  ('/local-data/DA-paper/reachr/2015_04_16_RR_F2_ab3_1_EA_1mM_2.mat');
haz_data = [23 24];
[stim,resp,ParadigmNames,paradigm] = MechanismAnalysis_PrepData(datapath,haz_data,stim,resp,ParadigmNames,paradigm,0);
alldata(2).stim = stim;
alldata(2).resp = resp;
alldata(2).ParadigmNames = ParadigmNames;
alldata(2).paradigm = paradigm;
MechanismAnalysis_PlotGain(stim,resp,ParadigmNames,paradigm,0);

PrettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% Here is another neuron, where, instead of steadily supplying light, we flicker the light @100Hz in an effort to get more sustained firing from the light:


datapath = '/local-data/DA-paper/reachr/2015_04_28_RR_F1_ab3_2_EA_400uM.mat';
haz_data = [14 16];
[stim,resp,ParadigmNames,paradigm] = MechanismAnalysis_PrepData(datapath,haz_data,[],[],{},[],0);
datapath = '/local-data/DA-paper/reachr/2015_04_28_RR_F1_ab3_2_EA_400uM_2.mat';

haz_data = 21;
[stim,resp,ParadigmNames,paradigm] = MechanismAnalysis_PrepData(datapath,haz_data,stim,resp,ParadigmNames,paradigm,0);
alldata(3).stim = stim;
alldata(3).resp = resp;
alldata(3).ParadigmNames = ParadigmNames;
alldata(3).paradigm = paradigm;
MechanismAnalysis_PlotGain(stim,resp,ParadigmNames,paradigm,0);

PrettyFig();

if being_published
	snapnow
	delete(gcf)
end

%% Summary
% Here we summarise the results from the three neurons show above. In the following figure, the fold change in gain is plotted for each neuron as a function of increasing light stimulation. Note that the odour stimulus doesn't change much. 

alldata(1).stim(:,7) = []; % this is an initial transient
alldata(1).resp(:,7) = [];
alldata(1).paradigm(:,7) = [];

MechanismAnalysis_PlotSummary(alldata);
PrettyFig();

if being_published
	snapnow
	delete(gcf)
end


%% Fast Gain Control
% Can we observe fast gain control in this dataset? Does it change with supplemental light stimulation? In the following figure, we analyse fast gain control in the first two neurons (the third was too noisy, and the stimulus wasn't reproducible enough). Brighter colours mean more light stimulation. It looks like the degree of fast gain control seems to be increasing with increasing light stimulation. 

c =parula(4);
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
for i = 1:length(alldata)-1
	clear ph
	ph(4) = subplot(1,length(alldata)-1,i); hold on
	for j = 1:length(alldata(i).ParadigmNames)
		% fit a filter
		resp = (alldata(i).resp(:,alldata(i).paradigm==j));
		stim = (alldata(i).stim(:,alldata(i).paradigm==j));
		if width(resp)>1
			resp  =mean2(resp);
			stim = mean2(stim);
		
			resp(1:1e4) = [];
			stim(1:1e4) = [];
			time = 1e-3*(1:length(stim));
			[K, ~, filtertime_full] = FindBestFilter(stim,resp,[],'regmax=1;','regmin=1;','filter_length=1999;','offset=500;');
			filtertime_full = filtertime_full*1e-3;
			filtertime = 1e-3*(-200:900);
			K = interp1(filtertime_full,K,filtertime);
			pred = convolve(time,stim,K,filtertime);
			%trivial scaling
			ff = fit(pred(~isnan(pred)),resp(~isnan(pred)),'poly1');
			pred = ff.p1*pred + ff.p2;


			history_lengths = logspace(log10(.1),1,30);
			example_history_length = history_lengths(10);
			[~,~,~,~,~,~,handles]= GainAnalysisWrapper2('response',resp,'prediction',pred,'stimulus',stim,'time',time,'engine',@GainAnalysis5,'history_lengths',history_lengths,'example_history_length',example_history_length,'ph',ph,'engine',@GainAnalysis5,'use_cache',true);

			% remove extra bullshit
			delete(handles.horz_line)
			delete(handles.vert_line)
			delete(handles.green_line)
			delete(handles.red_line)

			set(handles.red_dots,'MarkerFaceColor',c(j,:),'MarkerEdgeColor',c(j,:),'CData',c(j,:))
			set(handles.green_dots,'MarkerFaceColor',c(j,:),'MarkerEdgeColor',c(j,:),'CData',c(j,:))
		end

	end
	
end

PrettyFig();

if being_published
	snapnow
	delete(gcf)
end



%% Responses of ORNs to light flicker vs. odour flicker
% In this section, we present both a light flicker and a odour flicker to the same ORN and compare the filters we back out of the data. 

datapath = '/local-data/DA-paper/reachr/2015_04_17_RR_F3_ab3_1_EA_1mM.mat';
haz_data = 14;
[stim,resp,ParadigmNames,paradigm] = MechanismAnalysis_PrepData(datapath,haz_data,[],[],{},[],0);
haz_data = 18;
[stim,resp,ParadigmNames,paradigm] = MechanismAnalysis_PrepData(datapath,haz_data,stim,resp,ParadigmNames,paradigm,1);

c = parula(length(ParadigmNames)+1);
figure('outerposition',[0 0 1200 600],'PaperUnits','points','PaperSize',[1200 1000]); hold on
subplot(1,2,1), hold on
for i = 1:length(paradigm)
	[y,x] = hist(resp(1e4:end,i),50);
	plot(x,y,'Color',c(paradigm(i),:))
end
set(gca,'YScale','log')
xlabel('Response (Hz)')
ylabel('Count')

% show filters for each case
subplot(1,2,2), hold on
clear l
K = [];
for i = 1:width(stim)
	this_stim = stim(:,i);
	this_resp = resp(:,i);
	this_stim(1:1e4,:) = [];
	this_resp(1:1e4,:) = [];
	[this_K, ~, filtertime_full] = FindBestFilter(this_stim,this_resp,[],'regmax=1;','regmin=1;','filter_length=1999;','offset=500;');
	filtertime_full = filtertime_full*1e-3;
	filtertime = 1e-3*(-200:900);
	this_K = interp1(filtertime_full,this_K,filtertime);
	K = [K; this_K];
	this_K = this_K/max(this_K);
	l(paradigm(i)) = plot(filtertime,this_K,'Color',c(paradigm(i),:));
end
legend(l,ParadigmNames)
xlabel('Filter Lag (s)')
ylabel('Filter Amplitude (norm)')

PrettyFig();


if being_published
	snapnow
	delete(gcf)
end

%%
% Why are the filter shapes different? It's probably because we don't drive the ORNs strong enough with light. I would bet that if we drove the ORNs stronger, the filter shape will be the same. 


%%
% Can one filter be used to predict the response of the other? 

figure('outerposition',[0 0 1200 700],'PaperUnits','points','PaperSize',[1200 700]); hold on
subplot(2,2,1), hold on
title('Same Filter')
tA = 1e-3*(1:length(resp));
plot(tA,resp(:,1),'k')
ylabel('Response to odour (Hz)')
fp = convolve(tA,stim(:,1),K(1,:),filtertime);
% fix trival scaling
temp = fit(fp(1e4:end-200),resp(1e4:end-200,1),'poly1');
fp = fp*temp.p1 + temp.p2;
l = plot(tA,fp,'r');
r = rsquare(fp(1e4:end),resp(1e4:end,1));
legend(l,strcat('r^2=',oval(r)))
set(gca,'XLim',[20 40],'YLim',[0 100])

subplot(2,2,2), hold on
title('Different Filter')
plot(tA,resp(:,1),'k')
fp = convolve(tA,stim(:,1),K(4,:),filtertime);
% fix trival scaling
temp = fit(fp(1e4:end-200),resp(1e4:end-200,1),'poly1');
fp = fp*temp.p1 + temp.p2;
l = plot(tA,fp,'r');
r = rsquare(fp(1e4:end),resp(1e4:end,1));
legend(l,strcat('r^2=',oval(r)))
set(gca,'XLim',[20 40],'YLim',[0 100])

subplot(2,2,3), hold on
ylabel('Response to light (Hz)')
plot(tA,resp(:,4),'k')
fp = convolve(tA,stim(:,4),K(4,:),filtertime);
% fix trival scaling
temp = fit(fp(1e4:end-200),resp(1e4:end-200,4),'poly1');
fp = fp*temp.p1 + temp.p2;
l = plot(tA,fp,'r');
r = rsquare(fp(1e4:end),resp(1e4:end,4));
legend(l,strcat('r^2=',oval(r)))
set(gca,'XLim',[20 40],'YLim',[0 100])

subplot(2,2,4), hold on
plot(tA,resp(:,4),'k')
fp = convolve(tA,stim(:,4),K(1,:),filtertime);
% fix trival scaling
temp = fit(fp(1e4:end-200),resp(1e4:end-200,4),'poly1');
fp = fp*temp.p1 + temp.p2;
l = plot(tA,fp,'r');
r = rsquare(fp(1e4:end),resp(1e4:end,4));
legend(l,strcat('r^2=',oval(r)))
set(gca,'XLim',[20 40],'YLim',[0 100])

PrettyFig();

if being_published
	snapnow
	delete(gcf)
end


%%
% We know that ORNs adapt to steps of odour imperfectly, and this adaptation has two timescales. What about steps of light? Is the adaptation perfect? How many timescales are there? 

datapath = ('/local-data/DA-paper/reachr/2015_04_16_RR_F5_ab3_1_EA_1mM_2days.mat');
haz_data = 23;
[stim,resp,ParadigmNames,paradigm] = MechanismAnalysis_PrepData(datapath,haz_data,[],[],{},[],1);
resp(55e3:end)=[];
tA = 1e-3*(1:length(resp));

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(tA,resp,'k')
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')

PrettyFig();

if being_published
	snapnow
	delete(gcf)
end


%%
% While it's hard to say if there are two timescales or one, what is clear is that the degree of adaptation (the ratio of the peak firing rate to the adapted firing rate) is much higher than with an odour step that generates an equivalent peak firing rate. (This is partly remedied by flashing light @ 100Hz). 

%% Responses of ORNs to light flicker with an odour background
% In the following section, we do the corollary of the experiment we did before. Here, we present a flickering light stimulus, and present an odour background on top. 
clear alldata

datapath =  ('/local-data/DA-paper/reachr/2015_04_17_RR_F2_ab3_1_EA_1mM_4days.mat');
haz_data = [18:20];
[stim,resp,ParadigmNames,paradigm] = MechanismAnalysis_PrepData(datapath,haz_data,[],[],{},[],1);
alldata(1).stim = stim;
alldata(1).resp = resp;
alldata(1).ParadigmNames = ParadigmNames;
alldata(1).paradigm = paradigm;
MechanismAnalysis_PlotGain(stim,resp,ParadigmNames,paradigm,1);

PrettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% Here is another neuron:

datapath = ('/local-data/DA-paper/reachr/2015_04_17_RR_F3_ab3_4_EA_1mM.mat');
haz_data = [18:21];
[stim,resp,ParadigmNames,paradigm] = MechanismAnalysis_PrepData(datapath,haz_data,[],[],{},[],1);
stim(55e3:end,:) = [];
resp(55e3:end,:) = [];
alldata(2).stim = stim;
alldata(2).resp = resp;
alldata(2).ParadigmNames = ParadigmNames;
alldata(2).paradigm = paradigm;
MechanismAnalysis_PlotGain(stim,resp,ParadigmNames,paradigm,1);

PrettyFig();

if being_published
	snapnow
	delete(gcf)
end

%% 
% Here is yet another neuron, that we acquired with the LFP (not shown):
datapath = ('/local-data/DA-paper/reachr/2015_05_18_RR_F2_ab3_2_EA.mat');
haz_data = [25 26 27 28];
[stim,resp,ParadigmNames,paradigm] = MechanismAnalysis_PrepData(datapath,haz_data,[],[],{},[],1);
stim(1:1e4,:) = [];
resp(1:1e4,:) = [];
alldata(3).stim = stim;
alldata(3).resp = resp;
alldata(3).ParadigmNames = ParadigmNames;
alldata(3).paradigm = paradigm;
MechanismAnalysis_PlotGain(stim,resp,ParadigmNames,paradigm,1);

PrettyFig();
if being_published
	snapnow
	delete(gcf)
end



%%
% Here is a summary of all the data:

MechanismAnalysis_PlotSummary(alldata);
PrettyFig();

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

