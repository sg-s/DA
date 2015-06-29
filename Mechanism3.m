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


%% Responses to Odour Flicker with Light Backgrounds
% In the following figures, we present flickering odour stimuli to and record their responses to the odour with and without an additional light activation. 

allfiles = dir('/local-data/DA-paper/reachr/odour-flicker-light-background/2015_*.mat');
odour_flicker = [];
for i = 1:length(allfiles)
	datapath  = strcat('/local-data/DA-paper/reachr/odour-flicker-light-background/',allfiles(i).name);
	[stim,resp,ParadigmNames,paradigm] = MechanismAnalysis_PrepData(datapath,[],[],[],{},[],0);

	odour_flicker(i).stim = stim;
	odour_flicker(i).resp = resp;
	odour_flicker(i).ParadigmNames = ParadigmNames;
	odour_flicker(i).paradigm = paradigm;
	

	
end

%%
% In the following figure, we show that gain decreases as we increase the background light stimulation. The subplot on the left shows that the odour flickering stimulus that we deliver to the neuron is the same, no matter which paradigm we test, and that the odour stimulus doesn't change when we add the light. In the middle plot, we plot the firing rate trajectory without light background on the x-axis, and plot the firing rate trajectory with light background on the y-axis. The dotted line indicates unity. If the gain of the neuron changes with light background stimulation, the cloud of points will have a slope different from 1. If the gain doesn't change, but the mean firing rate is different, the cloud of points will move vertically up or down in the plot. In the third plot, we plot the mean firing rate in each paradigm (normalised to no background) on the y-axis, and the relative gain on the x-axis. 

MechanismAnalysis_PlotGain(odour_flicker(1).stim,odour_flicker(1).resp,odour_flicker(1).ParadigmNames,odour_flicker(1).paradigm,0);
PrettyFig();
if being_published
	snapnow
	delete(gcf)
end

%%
% Here is another neuron, just as before:

MechanismAnalysis_PlotGain(odour_flicker(4).stim,odour_flicker(4).resp,odour_flicker(4).ParadigmNames,odour_flicker(4).paradigm,0);
PrettyFig();
if being_published
	snapnow
	delete(gcf)
end



%% Responses of ORNs to light flicker with an odour background
% In the following section, we do the corollary of the experiment we did before. Here, we present a flickering light stimulus, and present an odour background on top. Each of the following figures shows the response from one ab3A neuron in a w; 22a-GAL4/+; UAS-ReaChR/+ fly. 
light_flicker = [];
allfiles = dir('/local-data/DA-paper/reachr/light-flicker-odour-background/2015_*.mat');
for i =1:length(allfiles)
	datapath  = strcat('/local-data/DA-paper/reachr/light-flicker-odour-background/',allfiles(i).name);
	[stim,resp,ParadigmNames,paradigm] = MechanismAnalysis_PrepData(datapath,[],[],[],{},[],1);

	light_flicker(i).stim = stim;
	light_flicker(i).resp = resp;
	light_flicker(i).ParadigmNames = ParadigmNames;
	light_flicker(i).paradigm = paradigm;

end

%%
% Here is an example neuron:
MechanismAnalysis_PlotGain(light_flicker(2).stim,light_flicker(2).resp,light_flicker(2).ParadigmNames,light_flicker(2).paradigm,1);
PrettyFig();
if being_published
	snapnow
	delete(gcf)
end

%%
% Here is another neuron, just as before:

MechanismAnalysis_PlotGain(light_flicker(4).stim,light_flicker(4).resp,light_flicker(4).ParadigmNames,light_flicker(4).paradigm,1);
PrettyFig();
if being_published
	snapnow
	delete(gcf)
end

%%
% We now summarise all the data, and compare with the odour flicker: 


figure('outerposition',[0 0 700 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
GainPhasePlot(odour_flicker,'b',1);
GainPhasePlot(light_flicker,'r',1);
xlabel('Relative Gain')
ylabel('Mean Firing rate (norm)')
PrettyFig();
if being_published
	snapnow
	delete(gcf)
end

%%
% This looks very unclear. Why is there so much variability? Here, we plot the response histograms of each point to see if there are some pattern which can clarify this picture:

figure('outerposition',[0 0 1400 800],'PaperUnits','points','PaperSize',[1400 800]); hold on
ncols = max([length(light_flicker) length(odour_flicker)]);
c = parula(5);
for i = 1:length(odour_flicker)
	subplot(2,ncols,i), hold on
	for j = 1:length(unique(odour_flicker(i).paradigm))
		this_resp = mean2(odour_flicker(i).resp(:,odour_flicker(i).paradigm == j));
		if isscalar(this_resp)
			this_resp = (odour_flicker(i).resp(:,odour_flicker(i).paradigm == j));
		end
		this_resp(isnan(this_resp)) = [];
		[y,x] = hist(this_resp,50);
		y = y/length(y); y = y/50;
		plot(x,y,'Color',c(j,:));
		
	end
end

for i = 1:length(light_flicker)
	subplot(2,ncols,i+ncols), hold on
	for j = 1:length(unique(light_flicker(i).paradigm))
		this_resp = mean2(light_flicker(i).resp(:,light_flicker(i).paradigm == j));
		if isscalar(this_resp)
			this_resp = (light_flicker(i).resp(:,light_flicker(i).paradigm == j));
		end
		this_resp(isnan(this_resp)) = [];
		[y,x] = hist(this_resp,50);
		y = y/length(y); y = y/50;
		plot(x,y,'Color',c(j,:));
		
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

