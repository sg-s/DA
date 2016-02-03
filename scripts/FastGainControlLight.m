% 
% 
% created by Srinivas Gorur-Shandilya at 1:53 , 02 February 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

%% Fast Gain Control w/o receptors?
% In this document, we try to see if we can see evidence of fast gain control in our data where we stimulate ORNs with light, but no odour. This document uses the new ORNData class, and all plots are generated using overloaded plot functions defined within this class. 

pHeader;

% load the data
if ~exist('light_data','var')
	load('/Users/sigbhu/code/da/data/Chrimson_data.mat') 
	% specify that we only want to use data from 5 to 55s
	temp = zeros(length(light_data.stimulus),1);
	temp(5e3:55e3) = 1;
	light_data.use_this_segment = logical(temp);

	% back out filters trial-wise
	light_data = backOutFilters(light_data);

end

%%
% Among all the data we have where we stimulate ORNs with just light, we pick stimulus paradigms where stimulus variance is large, as this is most likely to have evidence for fast gain control. Also, since this data was originally collected for a different purpose (to see if gain changes with mean stimulus), we have lots of experiments with lots of different stimuli, and it is important to only analyse groups of data where the stimulus is the same. 

%%
% In the following figure, we look at four different things in the data. First, we plot the firing rate vs. the projected stimulus, and look at the $r^2$ to ensure that our linear models are doing a good job. Second, we plot the distribution of instantaneous gain to see how much the gain varies over this experiment. Having established that the gain does vary a bit, we then plot the inst. gain vs. the mean stimulus in some preceding history length, for three different history lengths (10ms, 500ms, 10s). Finally, we plot the Spearman correlation vs. the history length to see if there is any window in which the mean stimulus predicts the inst. gain well. 

do_this_paradigm = 2;
light_data.use_these_trials = (light_data.paradigm == do_this_paradigm);
light_data = computeInstGain(light_data);

figure('outerposition',[0 0 1300 400],'PaperUnits','points','PaperSize',[1300 400]); hold on
ax(1) = subplot(1,4,1); hold on
ax(2) = subplot(1,4,2); hold on
ax(3) = subplot(1,4,3); hold on
ax(4) = subplot(1,4,4); hold on

% show the prediction vs. the data
plot(light_data,ax(1),'ioCurve.firing_rate','data_bin_type','dots','plot_type','mean','showr2',true);

% show the distribution of inst. gain
plot(light_data,ax(2),'pdf.inst_gain_firing','nbins',100);

plot(light_data,ax(3:4),'instGainAnalysis.firing_rate.mu','history_lengths',logspace(-2,1,30)*1e3,'history_length',[10 500 1e4],'data_bin_type','dots');
set(ax(3),'XScale','log','YScale','log','XLim',[1e-5 1e2])
prettyFig('fs=13;','FixLogX=true;','lw=1.3;','plw=1.2;')

if being_published
	snapnow
	delete(gcf)
end

%
% We now repeat this analysis, but for another stimulus paradigm:

do_this_paradigm = 3;
light_data.use_these_trials = (light_data.paradigm == do_this_paradigm);
light_data = computeInstGain(light_data);

figure('outerposition',[0 0 1300 400],'PaperUnits','points','PaperSize',[1300 400]); hold on
ax(1) = subplot(1,4,1); hold on
ax(2) = subplot(1,4,2); hold on
ax(3) = subplot(1,4,3); hold on
ax(4) = subplot(1,4,4); hold on

% show the prediction vs. the data
plot(light_data,ax(1),'ioCurve.firing_rate','data_bin_type','dots','plot_type','mean','showr2',true);

% show the distribution of inst. gain
plot(light_data,ax(2),'pdf.inst_gain_firing','nbins',100);

plot(light_data,ax(3:4),'instGainAnalysis.firing_rate.mu','history_lengths',logspace(-2,1,30)*1e3,'history_length',[10 500 1e4],'data_bin_type','dots');
set(ax(3),'XScale','log','YScale','log','XLim',[1e-5 1e2])
prettyFig('fs=13;','FixLogX=true;','lw=1.3;','plw=1.2;')

if being_published
	snapnow
	delete(gcf)
end



%% Version Info
%
pFooter;