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

figure('outerposition',[0 0 1400 800],'PaperUnits','points','PaperSize',[1400 800]); hold on
for i = 1:6
	ax(i) = subplot(2,3,i); hold on
end
hl = [50 300 1e4];
for i = 3:5
	plot_handles = plot(light_data,[ax(i) ax(6)],'instGainAnalysis.firing_rate.mu','history_lengths',logspace(-2,1,30)*1e3,'history_length',hl(i-2),'nbins',300);
	title(ax(i),['\tau_{H} = ' oval(hl(i-2)) 'ms'])
	if i < 5
		delete(plot_handles(2).f2)
	end
	ylabel(ax(i),'Inst Gain (Hz/\muW)')
end
set(ax(3:5),'XScale','log','YScale','log','YLim',[.5 10])
set(ax(6),'YLim',[-1 0])

% show the prediction vs. the data
plot(light_data,ax(1),'ioCurve.firing_rate','data_bin_type','dots','plot_type','mean','showr2',true);

% show the distribution of inst. gain
plot(light_data,ax(2),'pdf.inst_gain_firing','nbins',100);

prettyFig('fs=14;','FixLogY=true;');


if being_published
	snapnow
	delete(gcf)
end

%% 
% We now repeat this analysis, but for another stimulus paradigm:

do_this_paradigm = 3;
light_data.use_these_trials = (light_data.paradigm == do_this_paradigm);
light_data = computeInstGain(light_data);

figure('outerposition',[0 0 1400 800],'PaperUnits','points','PaperSize',[1400 800]); hold on
for i = 1:6
	ax(i) = subplot(2,3,i); hold on
end
hl = [50 300 1e4];
for i = 3:5
	plot_handles = plot(light_data,[ax(i) ax(6)],'instGainAnalysis.firing_rate.mu','history_lengths',logspace(-2,1,30)*1e3,'history_length',hl(i-2),'nbins',300);
	title(ax(i),['\tau_{H} = ' oval(hl(i-2)) 'ms'])
	if i < 5
		delete(plot_handles(2).f2)
	end
	ylabel(ax(i),'Inst Gain (Hz/\muW)')
end
set(ax(3:5),'XScale','log','YScale','log','YLim',[.5 10])
set(ax(6),'YLim',[-1 0])

% show the prediction vs. the data
plot(light_data,ax(1),'ioCurve.firing_rate','data_bin_type','dots','plot_type','mean','showr2',true);

% show the distribution of inst. gain
plot(light_data,ax(2),'pdf.inst_gain_firing','nbins',100);

prettyFig('fs=14;','FixLogY=true;');

if being_published
	snapnow
	delete(gcf)
end

%% Gain control beyond a static nonlinearity? 
% In this section we fit a full LN model to the data and check if there is gain control beyond a static nonlinearity. 


do_this_paradigm = 2;
light_data.use_these_trials = (light_data.paradigm == do_this_paradigm);
light_data = computeInstGain(light_data,true);

figure('outerposition',[0 0 1400 800],'PaperUnits','points','PaperSize',[1400 800]); hold on
for i = 1:6
	ax(i) = subplot(2,3,i); hold on
end
hl = [50 300 1e4];
for i = 3:5
	plot_handles = plot(light_data,[ax(i) ax(6)],'instGainAnalysis.firing_rate.mu','history_lengths',logspace(-2,1,30)*1e3,'history_length',hl(i-2),'nbins',300);
	title(ax(i),['\tau_{H} = ' oval(hl(i-2)) 'ms'])
	if i < 5
		delete(plot_handles(2).f2)
	end
	ylabel(ax(i),'Inst Gain (Hz/\muW)')
end
set(ax(3:5),'XScale','log','YScale','log','YLim',[.5 10])
set(ax(6),'YLim',[-1 .5])

% show the prediction vs. the data
plot(light_data,ax(1),'ioCurve.firing_rate','data_bin_type','dots','plot_type','mean','showr2',true);

% show the distribution of inst. gain
plot(light_data,ax(2),'pdf.inst_gain_firing','nbins',100);

prettyFig('fs=14;','FixLogY=true;');

if being_published
	snapnow
	delete(gcf)
end

%%
% Now we repeat it for another dataset.

do_this_paradigm = 3;
light_data.use_these_trials = (light_data.paradigm == do_this_paradigm);
light_data = computeInstGain(light_data,true);

figure('outerposition',[0 0 1400 800],'PaperUnits','points','PaperSize',[1400 800]); hold on
for i = 1:6
	ax(i) = subplot(2,3,i); hold on
end
hl = [50 300 1e4];
for i = 3:5
	plot_handles = plot(light_data,[ax(i) ax(6)],'instGainAnalysis.firing_rate.mu','history_lengths',logspace(-2,1,30)*1e3,'history_length',hl(i-2),'nbins',300);
	title(ax(i),['\tau_{H} = ' oval(hl(i-2)) 'ms'])
	if i < 5
		delete(plot_handles(2).f2)
	end
	ylabel(ax(i),'Inst Gain (Hz/\muW)')
end
set(ax(3:5),'XScale','log','YScale','log','YLim',[.5 10])
set(ax(6),'YLim',[-1 .5])

% show the prediction vs. the data
plot(light_data,ax(1),'ioCurve.firing_rate','data_bin_type','dots','plot_type','mean','showr2',true);

% show the distribution of inst. gain
plot(light_data,ax(2),'pdf.inst_gain_firing','nbins',100);

prettyFig('fs=14;','FixLogY=true;');

if being_published
	snapnow
	delete(gcf)
end


%% Version Info
%
pFooter;