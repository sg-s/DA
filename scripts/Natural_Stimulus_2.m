% Natural_Stimulus_2.m
%
% created by Srinivas Gorur-Shandilya at 11:20 , 10 February 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

%% Naturalistic Stimulus, take 2
% Here we document attempts to re-do the naturalistic stimulus, this time while also recording the LFP. 

pHeader;

path_name = '/local-data/DA-paper/natural-flickering/with-lfp/ab3';
[PID, LFP, fA, paradigm,orn,fly] = consolidateData(path_name,1);

% remove crap trials
rm_this = max(fA) == 0 | isnan(sum(LFP));
PID(:,rm_this) = [];
fA(:,rm_this) = [];
LFP(:,rm_this) = [];
paradigm(:,rm_this) = [];
fly(rm_this) = [];
orn(rm_this) = [];

% remove baselines from LFP and PID
for i = 1:width(PID)
	temp = nanmean(PID(1:4e3,i));
	PID(:,i) = PID(:,i) - temp;
	temp = nanmean(LFP(1:4e3,i));
	LFP(:,i) = LFP(:,i) - temp;
end

% band pass all the LFP
filtered_LFP = LFP;
for i = 1:width(LFP)
	filtered_LFP(:,i) = filtered_LFP(:,i) - fastFiltFilt(ones(1e4,1),1e4,filtered_LFP(:,i));
	filtered_LFP(:,i) = filtered_LFP(:,i)*10; % to get the units right, now in mV
end


%% Stimulus
% In this section we plot the stimulus we used. Because we're using Alicat MFCs rather than Aalborgs, this stimulus is different from the one originally used, despite an identical configuration of MFCs, valves and control signals. As you can see, this signal is denser, but still seems "naturalistic", with broad variation in the stimulus. Also, since we're using the Alicats, the stimulus is extremely reproducible trial-to-trial. The following figure shows the average of all the data, with a closeup to see the trial-to-trial variability. 

time = 1e-3*(1:length(PID));

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,5,1:4), hold on
plot(time,nanmean(PID,2),'k')
xlabel('Time (s)')
ylabel('Stimulus (V)')
title('Mean Stimulus')

subplot(1,5,5), hold on
a = 14e3; z = 18e3; ss = 10;
plot(time(a:ss:z),PID(a:ss:z,:))
xlabel('Time (s)')
title([oval(width(PID)) ' trials'])

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%% LFP and Firing Data
% The following figure shows how the raw and filtered LFP, together with the firing rate responses, averaged across each ORN.

figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
subplot(3,1,1), hold on
for i = 1:max(orn)
	plot(time,10*nanmean(LFP(:,orn==i),2))
end
set(gca,'XLim',[0 70])
ylabel('Raw LFP (mV)')
subplot(3,1,2), hold on
for i = 1:max(orn)
	plot(time,nanmean(filtered_LFP(:,orn==i),2))
end
set(gca,'XLim',[0 70])

ylabel('Filtered LFP (mV)')
subplot(3,1,3), hold on
for i = 1:max(orn)
	plot(time,nanmean(fA(:,orn==i),2))
end
set(gca,'XLim',[0 70])

ylabel('ORN Response (Hz)')
xlabel('Time (s)')
suptitle('Neuron-wise LFP and Firing Rates')
prettyFig('fs=14;');

if being_published
	snapnow
	delete(gcf)
end

%% Filters and Nonlinearities 
% In this section we back out filters for the LFP and the firing rate on a neuron-by-neuron basis. The following figure show the filters and and the residuals for the PID -> LFP transformation, for each ORN. Note that the LFP filter is simply integrating, and the nonlinearity is a simple curve.

% first convert everything to ORNData class
if ~exist('od','var')
	load('/local-data/DA-paper/natural-flickering/with-lfp/ab3/ORNData.mat','od')
else
	% od = ORNData;
	% for i = 1:max(orn)
	% 	disp(i)
	% 	od(i).stimulus = PID(:,orn==i);
	% 	od(i).LFP = filtered_LFP(:,orn==i);
	% 	od(i).firing_rate = fA(:,orn==i);

	% 	temp = zeros(7e4,1);
	% 	temp(15e3:55e3) = true;
	% 	od(i).use_this_segment = temp;
	% 	od(i).regularisation_factor = 1;	% this will automatically get the filters
	% end
end


for i = 1:max(orn)
	disp(i)
	od(i) = computeInstGain(od(i));
end

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
for i = 1:length(od)
	ax(1) = subplot(2,length(od),i); hold on; title(['ORN ' oval(i)])
	ax(2) = subplot(2,length(od),i+length(od)); hold on
	plot(od(i),ax,'LN.LFP','data_bin_type','dots','plot_type','mean','nbins',10)
end

prettyFig('fs=14;');

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

%% Inst. Gain Control
% Now we analyse the inst. gain in the LFP and the firing rate. 

for i = 1:max(orn)
	% first do LFP
	figure('outerposition',[0 0 1000 800],'PaperUnits','points','PaperSize',[1000 800]); hold on
	subplot(2,3,1), hold on
	plot(od(i),gca,'pdf.inst_gain_LFP','nbins',100);

	clear ax
	ax(1) = subplot(2,3,2); hold on
	ax(2) = subplot(2,3,3); hold on
	plot(od(i),ax,'instGainAnalysis.LFP.mu','history_lengths',logspace(-2,1,30)*1e3,'data_bin_type','dots');
	set(ax(1),'XScale','log','YScale','log')

	subplot(2,3,4), hold on
	plot(od(i),gca,'pdf.inst_gain_firing','nbins',100);

	clear ax
	ax(1) = subplot(2,3,5); hold on
	ax(2) = subplot(2,3,6); hold on
	plot(od(i),ax,'instGainAnalysis.firing_rate.mu','history_lengths',logspace(-2,1,30)*1e3,'data_bin_type','dots');
	set(ax(1),'XScale','log','YScale','log')

	prettyFig('fs=14;','FixLogX=true;');

	if being_published
		snapnow
		delete(gcf)
	end
end

%%
% Now we repeat the inst. gain analysis, but accounting for the nonlinearity:

warning off
for i = 1:max(orn)
	disp(i)
	od(i) = computeInstGain(od(i),true);
end
warning on


for i = 1:max(orn)
	% first do LFP
	figure('outerposition',[0 0 1000 800],'PaperUnits','points','PaperSize',[1000 800]); hold on
	subplot(2,3,1), hold on
	plot(od(i),gca,'pdf.inst_gain_LFP','nbins',100);

	clear ax
	ax(1) = subplot(2,3,2); hold on
	ax(2) = subplot(2,3,3); hold on
	plot(od(i),ax,'instGainAnalysis.LFP.mu','history_lengths',logspace(-2,1,30)*1e3,'data_bin_type','dots');
	set(ax(1),'XScale','log','YScale','log')

	subplot(2,3,4), hold on
	plot(od(i),gca,'pdf.inst_gain_firing','nbins',100);

	clear ax
	ax(1) = subplot(2,3,5); hold on
	ax(2) = subplot(2,3,6); hold on
	plot(od(i),ax,'instGainAnalysis.firing_rate.mu','history_lengths',logspace(-2,1,30)*1e3,'data_bin_type','dots');
	set(ax(1),'XScale','log','YScale','log')

	prettyFig('fs=14;','FixLogX=true;');

	if being_published
		snapnow
		delete(gcf)
	end
end


pFooter;

