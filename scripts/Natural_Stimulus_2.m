% Natural_Stimulus_2.m
%
% created by Srinivas Gorur-Shandilya at 11:20 , 10 February 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

%% Naturalistic Stimulus, take 2
% Here we document attempts to re-do the naturalistic stimulus, this time while also recording the LFP. 

pHeader;

% load the data
od = ORNData;
od = readData(od,'/local-data/DA-paper/natural-flickering/with-lfp/ab3/');


%% Stimulus
% In this section we plot the stimulus we used. Because we're using Alicat MFCs rather than Aalborgs, this stimulus is different from the one originally used, despite an identical configuration of MFCs, valves and control signals. As you can see, this signal is denser, but still seems "naturalistic", with broad variation in the stimulus. Also, since we're using the Alicats, the stimulus is extremely reproducible trial-to-trial. The following figure shows the average of all the data, with a closeup to see the trial-to-trial variability. 

time = 1e-3*(1:length(od(1).stimulus));

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
plot(time,[od.stimulus])
xlabel('Time (s)')
ylabel('Stimulus (V)')
title([oval(width([od.stimulus])) ' trials'])
set(gca,'XLim',[0 70])
prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%% LFP and Firing Data
% The following figure shows the stimulus, the LFP and the firing rate responses, averaged across each ORN.

figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
subplot(3,1,1), hold on
c = parula(length(od)+1);
for i = 1:length(od)
	plot(time,nanmean(od(i).stimulus,2),'Color',c(i,:))
end
set(gca,'XLim',[0 70])
ylabel('Stimulus (V)')
subplot(3,1,2), hold on
ylabel('\DeltaLFP (mV)')
for i = 1:length(od)
	plot(time,nanmean(od(i).LFP,2),'Color',c(i,:))
end
set(gca,'XLim',[0 70])
subplot(3,1,3), hold on
for i = 1:length(od)
	plot(time,nanmean(od(i).firing_rate,2),'Color',c(i,:))
end
set(gca,'XLim',[0 70])
ylabel('ORN Response (Hz)')
xlabel('Time (s)')
suptitle('Neuron-wise LFP and Firing Rates')
prettyFig('fs=14;','plw=1.5;');

if being_published
	snapnow
	delete(gcf)
end

%% Filters and Nonlinearities 
% In this section, we back out filters for the LFP and the firing rate. 


figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
for i = 1:length(od)
	ax(1) = subplot(2,length(od),i); hold on; title(['ORN ' oval(i)])
	ax(2) = subplot(2,length(od),i+length(od)); hold on
	plot(od(i),ax,'LN.LFP','data_bin_type','dots','plot_type','mean','nbins',10)
end


if being_published
	snapnow
	delete(gcf)
end

%%
% Now we look at the LN models for the firing rate.

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
for i = 1:length(od)
	ax(1) = subplot(2,length(od),i); hold on; title(['ORN ' oval(i)])
	ax(2) = subplot(2,length(od),i+length(od)); hold on
	plot(od(i),ax,'LN.firing_rate','data_bin_type','dots','plot_type','mean','nbins',10)
end

prettyFig('fs=14;');

if being_published
	snapnow
	delete(gcf)
end

%%
% We see that there is a weird non-monotonicity in the plot. We can verify that this actually exists in the data by comparing the projected stimulus to the firing rate directly: 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
x = nanmean(od(6).firing_projected,2);
y = nanmean(od(6).firing_rate,2);
ax = plotyy(time,x,time,y);
xlabel('Time (s)')
set(ax(1),'YLim',[-0.1 2],'XLim',[10 20])
set(ax(2),'YLim',[0 150],'XLim',[10 20])
ylabel(ax(1),'Firing Rate (Hz)')
ylabel(ax(2),'Proj. Stimulus (V)')
box off
prettyFig('fs=14;');

if being_published
	snapnow
	delete(gcf)
end

%% Fitting DA Models to the Data
% In this section, we fit DA models to the data on a neuron-by-neuron basis, to determine the timescale of fast gain control. The motivation for this analysis is two fold: if there is fast gain control, then DA models should outperform LN models. Also, the timescale of the gian filter in the DA model should be the timescale of the fast gain control. 
p.   s0 = -0.0935;
p.  n_z = 1.1328;
p.tau_z = 29.5858;
p.  n_y = 4.4375;
p.tau_y = 10.0293;
p.    C = 0.0047;
p.    A = 201.3815;
p.    B = 4.1471;
for j = 1:10
	for i = 1:length(od)
		clear d
		d.response = nanmean(od(i).firing_rate,2);
		d.stimulus = nanmean(od(i).stimulus,2);
		d.response(1:1e3) = NaN;
		p(i) = fitModel2Data(@DAModelv2,d,'p0',p(i),'nsteps',500,'make_plot',false);
	end
end

%% Dynamic Gain Control in LFP and the firing rate
% Now we analyse the inst. gain in the LFP and the firing rate. Having learnt our lesson about dynamic gain control being something that the LN model can't account for, we compute inst. gain vs. the output of a best-fit LN model. 

for i = 1:max(orn)
	od(i) = computeInstGain(od(i),true);
end

figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
for i = 1:6
	ax(i) = subplot(2,3,i); hold on
end
c = parula(max(orn)+1);

for i = 2:max(orn) % skip the first because it only has one trial
	% first do LFP inst gain distribution
	plot_handles = plot(od(i),ax(1),'pdf.inst_gain_LFP','nbins',100,'plot_type','mean');
	set(plot_handles.line(1),'Color',c(i,:))

	% show inst. gain analysis on LFP
	plot_handles = plot(od(i),ax(2:3),'instGainAnalysis.LFP.mu','history_lengths',logspace(-2,1,30)*1e3,'nbins',19);
	plot_handles(1).lines.Color = c(i,:);
	plot_handles(2).f2.Color = c(i,:);

	% now do firing inst. gain distribution
	plot_handles = plot(od(i),ax(4),'pdf.inst_gain_firing','nbins',100,'plot_type','mean');
	set(plot_handles.line(1),'Color',c(i,:))

	% show inst. gain analysis on firing rate
	plot_handles = plot(od(i),ax(5:6),'instGainAnalysis.firing_rate.mu','history_lengths',logspace(-2,1,30)*1e3,'nbins',19);
	plot_handles(1).lines.Color = c(i,:);
	plot_handles(2).f2.Color = c(i,:);
end
set(ax(2),'XScale','log','YScale','log','YLim',[.9 11])
set(ax(5),'XScale','log','YScale','log','YLim',[.9 11])
ylabel(ax(2),'Inst. Gain (norm)')
ylabel(ax(5),'Inst. Gain (norm)')
xlabel(ax(1),'Inst. gain (norm)')
xlabel(ax(4),'Inst. gain (norm)')
set(ax([1 4]),'XTick',[1e-2 1e-1 1e0 1e1 1e2])
title(ax(1),'LFP')
title(ax(4),'Firing Rate')
prettyFig('fs=14;','FixLogX=true;');

if being_published
	snapnow
	delete(gcf)
end

%% Whiff-based analysis
% In this section, we attempt to look at the gain on a whiff-by-whiff basis


R = nanmean(od(2).firing_rate,2);
[ons,offs] = computeOnsOffs(R>10);
rm_this = (offs-ons)>500;
ons(rm_this) = [];
offs(rm_this) = [];

figure, hold on
for i = 1:length(ons)
	plot(fp(ons(i):offs(i)),R(ons(i):offs(i)),'k.')
	pause(.5)
end


% for orn # 2
    % s0: -0.0148
    %   n_z: 2.2344
    % tau_z: 101.2500
    %   n_y: 8.2812
    % tau_y: 4.5752
    %     C: 0.7019
    %     A: 844.8893
    %     B: 10.6797


pFooter;

