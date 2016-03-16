% fast_gain_control_LFP.m
% 
% created by Srinivas Gorur-Shandilya at 7:15 , 15 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.



pHeader;

%% Fast Gain Control in LFPs?
% In this document we see if we can detect dynamical gain control in the LFP signal. The analysis is exactly the same as we did for the firing rates.


% load the data
if ~exist('od','var')
	p = '/local-data/DA-paper/natural-flickering/with-lfp/ab3/';
	od = raw2ORNData(p,'filter_LFP',true);

	% specify interesting range
	uts = zeros(length(od(1).stimulus),1);
	uts(10e3:end-10e3) = true;
	for i = 1:length(od)
		od(i).use_this_segment = uts;
	end

	% back out all filters
	od = backOutFilters(od,'all');
end

figure('outerposition',[0 0 1500 800],'PaperUnits','points','PaperSize',[1500 800]); hold on
clear ax
ax(1) = subplot(2,2,1); hold on
ax(2) = subplot(2,2,2); hold on
for i = 3:6
	ax(i) = subplot(2,4,2+i); hold on
end

% show stimulus
t = (1:length(od(1).stimulus))*1e-3;
plot(ax(1),t,nanmean([od.stimulus],2),'k')
xlabel(ax(1),'Time (s)')
ylabel(ax(1),'Stimulus (V)')
set(ax(1),'XLim',[0 70])

% show response
for i = 1:length(od)
	plot(ax(2),t,nanmean([od(i).LFP],2),'Color',[0 0 0 0.1])
end
plot(ax(2),t,nanmean([od.LFP],2),'Color','k')
xlabel(ax(2),'Time (s)')
ylabel(ax(2),'LFP (mV)')
set(ax(2),'XLim',[0 70])

% do the analysis of fast gain control
example_data = od(2);
orn_data = ORNData;
orn_data.stimulus = nanmean(example_data.stimulus,2);
orn_data.LFP = nanmean(example_data.LFP,2);
stim_on = false(length(example_data.stimulus),1);
stim_on(10e3:end-5e3) = true; 
orn_data.use_this_segment = stim_on;

options = fastGainControlAnalysis;
options.data = 'LFP';
r2_plot_data = fastGainControlAnalysis(ax(3:end),orn_data,options);

legend(r2_plot_data.l,['Peak = ' oval(r2_plot_data.peak_tau) 'ms'])
prettyFig('fs=18;','FixLogX=true;')

if being_published	
	snapnow	
	delete(gcf)
end

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end


%% Version Info
%
pFooter;


