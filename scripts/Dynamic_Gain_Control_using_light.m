% Dynamic_Gain_Control_using_light.m
% 
% created by Srinivas Gorur-Shandilya at 3:53 , 02 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

% this loads and assembles the data
% x = [.75:0.05:1.1 1.2:.1:3.6 3.8 4 4.2 4.5];
% x = [0:0.1:0.7 x];
% y = [0 4 12.8 24.1 36.2 48.1 60 70 95 114 134 151 167 184 201 219 236 252 269 283 299 314 328 343 357 370 383 394 404 413 419 422 422 421 419 418 417];
% y = [0*(0:0.1:0.7) y ];
% light_power_fit = fit(x(:),y(:),'spline');

% path_name = '/local-data/DA-paper/Chrimson/nat-stim-mimic/';
% [od] = raw2ORNData(path_name,light_power_fit);

% % remove the LFP
% for i = 1:size(od,1)
% 	for j = 1:size(od,2)
% 		od(i,j).LFP = [];
% 	end
% end

% % back out filters
% od = backOutFilters(od);

load('/local-data/DA-paper/Chrimson/nat-stim-mimic/combined_data.ORNData','od','-mat')


pHeader;

%% Dynamic Gain Control using Light
% In this document we check if we can see signs of dynamic gain control using just light stimulation, and no odour. The stimulus we present is chosen to mimic the naturalistic odour stimulus, and looks like this: 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(od(6,6),gca,'timeseries.stimulus','plot_type','mean');
set(gca,'YLim',[0 100],'XLim',[0 60])
ylabel('Light Power (\muW)')

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%%
% We now perform a excursion based analysis. This is how it works: first, we fit a LN model to the data. Deviations from the LN model are presumed to be associated with gain changes in the ORN. To measure gain changes, we first find all excursions in the dataset. Excursions are defined by deviations in the response (here, firing rate) exceeding the floor by 10%. 

%%
% In the following figure, we show plots of response vs. LN model prediction for three ORNs across 4 different stimuli, slightly scaled from one another. In the next figure, we show the gain during each excursion as a function of the mean stimulus in the past 500ms, and also plot the Spearman correlation coefficient as a function of history length. 

% now fit pwlinear functions to the resp. vs. proj. stimulus
od2 = fitNL(od);

f1 = figure('outerposition',[0 0 1300 800],'PaperUnits','points','PaperSize',[1300 800]); hold on
clear ax1 ax2
for i = 1:6
	ax1(i) = subplot(2,3,i); hold on
end
f2 = figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ax2 = subplot(1,2,1); hold on
ax3 = subplot(1,2,2); hold on

do_these = [6 6; 6 3; 6 1; 5 1; 2 3; 2 2];
for i = 1:length(do_these)
	a = do_these(i,1); b = do_these(i,2);
	plot(od2(a,b),[ax1(i) ax2 ax3],'excGainAnalysis.firing_rate.mu');
end

set(ax2,'YScale','log','YLim',[.5 20])
c = lines(7);
h = ax2.Children;
for i = 1:length(h)
	set(h(i),'Color',c(i,:),'Marker','o','MarkerSize',6,'MarkerFaceColor',c(i,:))
end

h = ax3.Children;
for i = 1:length(h)
	set(h(i),'Color',c(i,:),'Marker','o','MarkerSize',6,'MarkerFaceColor',c(i,:))
end

figure(f1)
prettyFig('fs=14;')

figure(f2)
prettyFig('fs=14;')


if being_published	
	snapnow	
	delete(f1)
	delete(f2)
end

%%
% While we do see signs of dynamical gain control here, what is weird is that it does not go to zero as the history length goes to zero. Why is that? 

%%
% One idea is that stimulus correlations are such that they prevent us from effectively determining the time window of gain control. It could work like this: imagine a long pulse of stimulus, with brief probes towards the end, followed by a long silence, followed by more probes. The ORN responds with lower gain to the first set of probes than the second, because it adapted to the stimulus over the last X ms. However, even the stimulus over the last Y ms is correlated with this gain change, where Y << X. 

%%
% To test if this is true, I generated some synthetic data using a DA model fit to these ORN responses, and driven using the same stimulus as used in these experiments. I then repeated this analysis, to see if I see a similar effect. 

clear p
p.tau_y = 15.6250;
p.tau_z = 257.5000;
p.  n_z = 2;
p.  n_y = 3.7500;
p.    A = 32.7500;
p.    B = 2.0312;
p.    C = 0.0812;
p.   s0 = -3.8125;
dd = ORNData;
dd.stimulus = nanmean(od(6,6).stimulus,2);
dd.firing_rate = abs(DAModelv2(dd.stimulus,p));
dd = backOutFilters(dd);
dd = fitNL(dd);

figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[1000 800]); hold on
clear ax
for i = 1:4
	ax(i) = subplot(2,2,i); hold on
end
plot(dd,ax(1),'Filter.firing_rate');
plot(dd,ax(2:end),'excGainAnalysis.firing_rate.mu');

ylabel(ax(2),'DA model output (Hz)')
set(ax(3),'YScale','log','YLim',[0.1 10])

prettyFig('fs=14;')

if being_published	
	snapnow	
	delete(gcf)
end

%%
% We see that in this synthetic data, we still see very low Spearman correlation coefficents even as the history length tends to zero, and especially for values within the timescale of the filter in the LN model.

%% Directly Estimating Gain timescales
% Since we define gain control in terms of deviations from a LN model, and since we believe gain occurs divisively over some recent timescale in a stimulus-driven manner, we can directly fit an adaptive gain model to the output of a LN model and check if it can account for any of the variance unexplained by a LN model. 

do_these = [6 6; 6 3; 6 1; 5 1; 2 3; 2 2];
% clear p p0 d
% p0.A = 2; p0.B = 0.05; p0.tau = 500; p0.n = 1;
% for i = 1:length(do_these)
% 	a = do_these(i,1); b = do_these(i,2);
% 	d.stimulus = nanmean(od2(a,b).firing_projected,2);
% 	d.response = nanmean(od2(a,b).firing_rate,2);
% 	p(i) = fitModel2Data(@adaptiveGainModel,d,'p0',p0);
% end
% save('.cache/adaptive_gain_model_fits_light_flicker.mat','p')

load('.cache/adaptive_gain_model_fits_light_flicker.mat','p')

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ax(1) = subplot(1,2,1); hold on
ax(2) = subplot(1,2,2); hold on
for i = 1:length(do_these)
	a = do_these(i,1); b = do_these(i,2);
	d.stimulus = nanmean(od2(a,b).firing_projected,2);
	d.response = nanmean(od2(a,b).firing_rate,2);
	r2LN = rsquare(d.stimulus,d.response);
	fp = adaptiveGainModel(d.stimulus,p(i));
	r2AG = rsquare(d.response,fp);
	plot(ax(1),r2LN,r2AG,'k+')
	plot(ax(2),100*(-r2LN+r2AG),p(i).tau*p(i).n,'k+')
end
plot(ax(1),[0 1],[0 1],'k--')
set(ax(1),'XLim',[.7 1],'YLim',[.7 1])
xlabel(ax(1),'r^2_{LN}')
ylabel(ax(1),'r^2_{adaptive gain}')

xlabel(ax(2),'Excess variance explained (%)')
ylabel(ax(2),'Timescale of gain control (ms)')

prettyFig('fs=18;')

if being_published	
	snapnow	
	delete(gcf)
end


%% Version Info
%
pFooter;




