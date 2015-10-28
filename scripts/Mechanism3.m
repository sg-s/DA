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

%% Kinetics
% First we study the response delays w.r.t to the PID with and without supplemental light stimulation. We estimate response delays using cross correlations. 

[PID, LFP, fA, paradigm, orn, fly, AllControlParadigms, paradigm_hashes, sequence] = consolidateData('/local-data/DA-paper/reachr/odour-flicker-light-background/',true);

% recompute paradigms sensibly 
light_V = NaN*paradigm;
for i = 1:length(paradigm)
	p = paradigm(i);
	if isempty(strfind(AllControlParadigms(p).Name,'LightFlicker')) & ~isempty(strfind(AllControlParadigms(p).Name,'Odour')) & ~isempty(strfind(AllControlParadigms(p).Name,'Light'))
		light_V(i) = str2double(AllControlParadigms(p).Name(strfind(AllControlParadigms(p).Name,'+')+1:strfind(AllControlParadigms(p).Name,'V')-1));
	end
end

% trash some trash
rm_this = false(width(PID),1);
for i = 1:width(PID)
	if any(isnan(PID(:,i)))
		rm_this(i) = true;
	end
	if max(PID(30e3:50e3,i)) - min(PID(30e3:50e3,i)) < .4
		rm_this(i) = true;
	end
end
PID(:,rm_this) = [];
fA(:,rm_this) = [];
paradigm(rm_this) = [];
light_V(rm_this) = [];

% remove some junk
light_V(light_V == 4) = NaN;
light_V(light_V == .8) = NaN;
light_V(light_V == 1.1) = NaN;
light_V(light_V == 1.2) = NaN;

% convert light voltage into power
load('/local-data/DA-paper/LED_calibration_data.mat','V','lux');
ff = fit(V(:),lux(:),'smoothingspline');
light_power = ff(light_V);
light_power(light_power<0) = 0;

% extract filters
a = 30e3; z = 55e3;
[K, prediction, gain, gain_err] = extractFilters(PID,fA,'a',a,'z',z);

% compute response delays (PID->firing)
tau = NaN(width(PID),1);
tau_LFP = NaN*tau;
for i = 1:width(PID)
	if ~any(isnan(K(:,i))) & mean(fA(a:z,i)) > 10
		x = PID(a:z,i);
		x = x - mean(x); x = x/std(x);
		y = fA(a:z,i);
		y = y - mean(y); y = y/std(y);
		temp = xcorr(x,y);
		[~,loc] = max(temp);
		tau(i) = z-a - loc;
	end
end 



figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
plot(light_power,tau,'k+')
ylabel('Response Delay (ms)')
xlabel('Supplemental Light (\muW)')

subplot(1,2,2), hold on
l = plot(gain,tau,'k+');
legend(l,['\rho= ' oval(spear(gain,tau))])
xlabel('Gain (Hz/V)')
ylabel('Response Delay (ms)')
prettyFig();

if being_published
	snapnow
	delete(gcf)
end

return
%% Gain Changes

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
all_light_power = nonnans(unique(light_power));
for i = 1:length(all_light_power)
	plot_this= find(light_power == all_light_power(i));
	errorbar(all_light_power(i),mean(gain(plot_this)),sem(gain(plot_this)),'k')
end

ylabel('Gain (Hz/V)')
xlabel('Supplemental Light (\muW)')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end


allfiles = dir('/local-data/DA-paper/reachr/odour-flicker-light-background/2015_*.mat');
odour_flicker = [];
for i = 1:length(allfiles)
	datapath  = strcat('/local-data/DA-paper/reachr/odour-flicker-light-background/',allfiles(i).name);
	[stim,resp,ParadigmNames,paradigm] = mechanismAnalysisPrepData(datapath,[],[],[],{},[],0);

	odour_flicker(i).stim = stim;
	odour_flicker(i).resp = resp;
	odour_flicker(i).ParadigmNames = ParadigmNames;
	odour_flicker(i).paradigm = paradigm;
end

%%
% In the following figure, we show that gain decreases as we increase the background light stimulation. The subplot on the left shows that the odour flickering stimulus that we deliver to the neuron is the same, no matter which paradigm we test, and that the odour stimulus doesn't change when we add the light. In the middle plot, we plot the firing rate trajectory without light background on the x-axis, and plot the firing rate trajectory with light background on the y-axis. The dotted line indicates unity. If the gain of the neuron changes with light background stimulation, the cloud of points will have a slope different from 1. If the gain doesn't change, but the mean firing rate is different, the cloud of points will move vertically up or down in the plot. In the third plot, we plot the mean firing rate in each paradigm (normalised to no background) on the y-axis, and the relative gain on the x-axis. 

mechanismAnalysisPlotGain(odour_flicker(1).stim,odour_flicker(1).resp,odour_flicker(1).ParadigmNames,odour_flicker(1).paradigm,0);
PrettyFig();
if being_published
	snapnow
	delete(gcf)
end

%%
% Here is another neuron, just as before:

mechanismAnalysisPlotGain(odour_flicker(4).stim,odour_flicker(4).resp,odour_flicker(4).ParadigmNames,odour_flicker(4).paradigm,0);
PrettyFig();
if being_published
	snapnow
	delete(gcf)
end

return

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
i=1;
MechanismAnalysis_PlotGain(light_flicker(i).stim,light_flicker(i).resp,light_flicker(i).ParadigmNames,light_flicker(i).paradigm,1);
PrettyFig();
if being_published
	snapnow
	delete(gcf)
end

%%
% We now summarise all the data, and compare with the odour flicker: 


figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
GainPhasePlot(odour_flicker,'b',0);
xlabel('Relative Gain')
ylabel('Mean Firing rate (Hz)')
title('Odour flicker, light background')
set(gca,'XLim',[0.3 1.3],'YLim',[0 65])
subplot(1,2,2), hold on
GainPhasePlot(light_flicker,'r',0);
xlabel('Relative Gain')
title('Light flicker, odour background')
set(gca,'XLim',[0.3 1.3],'YLim',[0 65])
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

