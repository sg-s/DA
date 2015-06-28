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

allfiles = dir('/local-data/DA-paper/reachr/odour-flicker-light-background/2015_*.mat');
% allfiles([2 3 4 10]) = [];
odour_flicker = [];
for i = 1:length(allfiles)
	datapath  = strcat('/local-data/DA-paper/reachr/odour-flicker-light-background/',allfiles(i).name);
	[stim,resp,ParadigmNames,paradigm] = MechanismAnalysis_PrepData(datapath,[],[],[],{},[],0);

	odour_flicker(i).stim = stim;
	odour_flicker(i).resp = resp;
	odour_flicker(i).ParadigmNames = ParadigmNames;
	odour_flicker(i).paradigm = paradigm;
	% MechanismAnalysis_PlotGain(stim,resp,ParadigmNames,paradigm,0);
	% PrettyFig();

	% if being_published
	% 	snapnow
	% 	delete(gcf)
	% end
end




%% Responses of ORNs to light flicker with an odour background
% In the following section, we do the corollary of the experiment we did before. Here, we present a flickering light stimulus, and present an odour background on top. Each of the following figures shows the response from one ab3A neuron in a w; 22a-GAL4/+; UAS-ReaChR/+ fly. 
allfiles = dir('/local-data/DA-paper/reachr/light-flicker-odour-background/2015_*.mat');
for i =1:length(allfiles)
	datapath  = strcat('/local-data/DA-paper/reachr/light-flicker-odour-background/',allfiles(i).name);
	[stim,resp,ParadigmNames,paradigm] = MechanismAnalysis_PrepData(datapath,[],[],[],{},[],1);

	alldata(i).stim = stim;
	alldata(i).resp = resp;
	alldata(i).ParadigmNames = ParadigmNames;
	alldata(i).paradigm = paradigm;
	% MechanismAnalysis_PlotGain(stim,resp,ParadigmNames,paradigm,1);
	% PrettyFig();

	% if being_published
	% 	snapnow
	% 	delete(gcf)
	% end

end

%%
% We now summarise all the data, and compare with the odour flicker: 


datapath =  ('/local-data/DA-paper/reachr/light-flicker-odour-background/2015_04_17_RR_F2_ab3_1_EA_1mM_4days.mat');
haz_data = [18:20];
[stim,resp,ParadigmNames,paradigm] = MechanismAnalysis_PrepData(datapath,haz_data,[],[],{},[],1);



%%
% Here is the 2nd neuron:

datapath = ('/local-data/DA-paper/reachr/light flicker + odour background/2015_04_17_RR_F3_ab3_4_EA_1mM.mat');
haz_data = [18:21];
[stim,resp,ParadigmNames,paradigm] = MechanismAnalysis_PrepData(datapath,haz_data,[],[],{},[],1);
stim(55e3:end,:) = []; % because we lost some data here
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
% Here is a 3rd neuron, that we acquired with the LFP (not shown):

datapath = ('/local-data/DA-paper/reachr/light flicker + odour background/2015_05_18_RR_F2_ab3_2_EA.mat');
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
% In this figure, we compare the slope and the offset of each data set, in the odour flicker and the light flicker case: 

figure('outerposition',[0 0 700 700],'PaperUnits','points','PaperSize',[700 700]); hold on
GainPhasePlot(odour_flicker,'b');
GainPhasePlot(alldata,'r')
ylabel('Mean Firing Rate (norm)')
xlabel('Gain (norm)')

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

