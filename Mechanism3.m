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

MechanismAnalysis_PlotGain(stim,resp,ParadigmNames,paradigm,0);

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


figure('outerposition',[0 0 600 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
% show filters for each case
c = parula(length(ParadigmNames)+1);
clear l
for i = 1:width(stim)
	this_stim = stim(:,i);
	this_resp = resp(:,i);
	this_stim(1:1e4,:) = [];
	this_resp(1:1e4,:) = [];
	[this_K, ~, filtertime_full] = FindBestFilter(this_stim,this_resp,[],'regmax=1;','regmin=1;','filter_length=1999;','offset=500;');
	filtertime_full = filtertime_full*1e-3;
	filtertime = 1e-3*(-200:900);
	this_K = interp1(filtertime_full,this_K,filtertime);
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
% Why are the filter shapes different? What does it mean?

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

% [~,loc]=max(resp);
% t = tA(loc:end);
% mint = min(t);
% t = t -min(t);
% ft=fittype('a*exp(-x./b1)+ a*exp(-x./b2) +c');
% fo = fitoptions(ft);
% fo.Robust = 'on';
% tf = mean(resp(end-100:end));
% pf = max(resp);  pf = pf - tf; pf = pf/2;
% fo.Lower = [pf-1e-3 1e-2 1e-2 tf-1e-3];
% fo.Upper = [pf+1e-3 2 20 tf+1e-3];
% fo.StartPoint = [pf .1 1 tf];
% ff = fit(t(:),resp(loc:end),ft,fo);
% plot(t+mint,ff(t),'r')

PrettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% While it's hard to say if there are two timescales or one, what is clear is that the degree of adaptation (the ratio of the peak firing rate to the adapted firing rate) is much higher than with an odour step that generates an equivalent peak firing rate. 

%% Responses of ORNs to light flicker with an odour background
% In the following section, we do the corollary of the experiment we did before. Here, we present a flickering light stimulus, and present an odour background on top. 

datapath =  ('/local-data/DA-paper/reachr/2015_04_17_RR_F2_ab3_1_EA_1mM_4days.mat');
haz_data = [18:20];
[stim,resp,ParadigmNames,paradigm] = MechanismAnalysis_PrepData(datapath,haz_data,[],[],{},[],1);
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
stim(55e4:end,:) = [];
resp(55e4:end,:) = [];
MechanismAnalysis_PlotGain(stim,resp,ParadigmNames,paradigm,1);

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

