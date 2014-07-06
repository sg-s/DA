% Mahmut_Data_Analyis.m
% 
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

% load data
load('/local-data/DA-paper/mahmut_data.mat')

%% 
% How do ORNs respond to non-Gaussian inputs? Real odor stimuli are characterised by long tails and non-gaussian statisitcs, with large whiffs of odor that occur in periods of relatively low signal. Such a stimulus has been generated here, and the responses of ORNs to these signals is analysed in this figure.

redo_bootstrap = 0;

td = 1;

%% Stimulus Characteristics 
% The following figure shows what the stimulus and the neuron's looks like. 

% first, remove the baseline from the PID
data(td).PID = data(td).PID - mean(data(td).PID(1:3000));
data(td).PID(data(td).PID < 0) = 0;

figure('outerposition',[0 0 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
subplot(2,1,1), hold on
plot(data(td).time,data(td).PID,'k');
ylabel('PID (V)')
title('ab3A response to ethyl acetate')
set(gca,'YLim',[-0.1 3.1])

subplot(2,1,2), hold on
plot(data(td).time,data(td).ORN,'k');
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')

PrettyFig;

%% Stimulus and Response Statistics 
% The following figure describes the statistics of the stimulus and the response. Left panel: Histograms of stimulus and response. Middle panel: Autocorrelation functions of the stimulus and the response. Right: Linear filter extracted from this dataset. 

% build a simple linear model
[K,~,filtertime] = FindBestFilter(data(td).PID(500:end),data(td).ORN(500:end),[],'filter_length=201;','min_cutoff = 0;');
data(td).K = K;
data(td).filtertime = filtertime*mean(diff(data(td).time));
data(td).LinearFit = mean(data(td).ORN) + convolve(data(td).time,data(td).PID,data(td).K,data(td).filtertime);
data(td).LinearFit(data(td).LinearFit < 0) = 0;

clear ph
figure('outerposition',[0 0 1300 400],'PaperUnits','points','PaperSize',[1300 400]); hold on
ph(1)=subplot(1,3,1); hold on
ph(2)=subplot(1,3,2); hold on
[act] = PlotDataStatistics(data,td,ph);

subplot(1,3,3), hold on;
plot(data(td).filtertime,data(td).K,'k','LineWidth',2)
set(gca,'XLim',[min(data(td).filtertime) max(data(td).filtertime)],'box','on')
xlabel('Filter Lag (s)')
title('Filter')
PrettyFig;

snapnow;
delete(gcf);

%% Variable Response
% The most intereesting thing about this dataset is that sometimes, for some whiffs, the ORN responds, but for others, it does not. The following figure shows a close-up of the data to illustrate this point. On the left, the traces show that ORNs respond only to some pulses, but not to others. It's not clear why. On the right, the plots show that ORN response to two almost identical pulses is very different, with the response to the second pulse much smaller than the response to the first. 


figure('outerposition',[0 0 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
subplot(2,2,1), hold on
plot(data(td).time,data(td).PID,'k');
ylabel('PID (V)')
set(gca,'YLim',[-0.1 3.1],'XLim',[25 35])

subplot(2,2,2), hold on
plot(data(td).time,data(td).PID,'k');
ylabel('PID (V)')
set(gca,'YLim',[-0.1 3.1],'XLim',[38 48])

subplot(2,2,3), hold on
plot(data(td).time,data(td).ORN,'k');
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
set(gca,'XLim',[25 35])

subplot(2,2,4), hold on
plot(data(td).time,data(td).ORN,'k');
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
set(gca,'XLim',[38 48])
PrettyFig;

%%
% Why is there so much variability in the response? The following figure shows a scatter plot of the stimulus and the response. 
figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(data(td).PID,data(td).ORN,'.k','MarkerSize',24)
ylabel('Firing Rate (Hz)')
xlabel('Stimulus (V)')
PrettyFig;

%%
% There is surprisingly little correlation between the input and the output. In particular, trajectories seem to cover all the space, instead of being confined to one curve. The r-square between the stimulus and the response is:

disp(rsquare(data(td).PID,data(td).ORN))

%%
% Let's attempt to find all the times when the stimulus is high, i.e., ten standard deviations above the baseline. This threshold ensures that we pick up all the whiffs, but ignore the blanks in between, as shown in the figure below. 

figure('outerposition',[0 0 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
plot(data(td).time, data(td).PID,'k'), hold on

hs = data(td).PID > 10*std(data(td).PID(1:3000));
plot(data(td).time(hs), data(td).PID(hs),'.r','MarkerSize',24), hold on
set(gca,'XLim',[28 36],'YLim',[-0.1 3.1])
xlabel('Time (s)')
ylabel('Stimulus (V)')
PrettyFig;

%%
% Now we break up the trace so that we can perform a whiff-by-whiff analysis of the response. The following figure shows the relationship between response and stimulus on a whiff-by-whiff basis. 

[whiff_ons,whiff_offs] = ComputeOnsOffs(hs);

whiff_durations = whiff_offs - whiff_ons;
whiff_ons(whiff_durations <20) = [];
whiff_offs(whiff_durations <20) = [];
whiff_durations = whiff_offs - whiff_ons;

whiff_stim_max  = NaN(1,length(whiff_ons));
whiff_resp_max  = NaN(1,length(whiff_ons));
whiff_stim_sum  = NaN(1,length(whiff_ons));
whiff_resp_sum  = NaN(1,length(whiff_ons));

for i = 1:length(whiff_ons)
	whiff_stim_sum(i) = sum(data(td).PID(whiff_ons(i):whiff_offs(i)));
	whiff_stim_max(i) = max(data(td).PID(whiff_ons(i):whiff_offs(i)));
	whiff_resp_sum(i) = sum(data(td).ORN(whiff_ons(i):whiff_offs(i)))*mean(diff(data(td).time));
	whiff_resp_max(i) = max(data(td).ORN(whiff_ons(i):whiff_offs(i)));
end
clear i

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
plot(whiff_stim_sum, whiff_resp_sum,'.r','MarkerSize',24), hold on
xlabel('Sum Whiff Stimulus (V s)')
ylabel('Sum ORN response (spikes)')

subplot(1,2,2), hold on
plot(whiff_stim_max, whiff_resp_max,'.r','MarkerSize',24), hold on
xlabel('Max Whiff Stimulus (V)')
ylabel('Max ORN response (Hz)')

PrettyFig;

%%
% The sum total of the odor delivered/whiff seems to correlate well with the number of spikes elicited in that whiff. The r-square of this is:

disp(rsquare(whiff_stim_sum,whiff_resp_sum));

%%
% The examination of the raw data strongly suggests that the response of the ORN seems to depend on the previous stimulus. (e.g., response to pulse at $t=44s$ is half that of a very similar pulse at $t=39s$)

%% DA Model fit to Data
% Can we fit a DA Model to this data? Does it explain the observed variability? The following figure shows the ORN firing rates and the best-fit DA Model. 

data(1).DAFit=DA_integrate2(data(1).PID,data(1).DAFitParam);

figure('outerposition',[0 0 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
subplot(2,2,1), hold on
plot(data(td).time,data(td).PID,'k');
ylabel('PID (V)')
set(gca,'YLim',[-0.1 3.1],'XLim',[25 35])

subplot(2,2,2), hold on
plot(data(td).time,data(td).PID,'k');
ylabel('PID (V)')
set(gca,'YLim',[-0.1 3.1],'XLim',[38 48])

subplot(2,2,3), hold on
plot(data(td).time,data(td).ORN,'k');
plot(data(1).time,data(1).DAFit,'r')
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
set(gca,'XLim',[25 35])

subplot(2,2,4), hold on
plot(data(td).time,data(td).ORN,'k');
plot(data(1).time,data(1).DAFit,'r')
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
set(gca,'XLim',[38 48])
PrettyFig;
legend Data DAFit

%%
% The r-square of the fit is:

disp(rsquare(data(1).ORN,data(1).DAFit))

return
	


	

	


	% do gain analysis
	figure('outerposition',[0 0 1000 400],'PaperUnits','points','PaperSize',[1300 400]); hold on
	clear ph
	ph(3)=subplot(1,2,1); hold on; 	axis square
	ph(4)=subplot(1,2,2); hold on;	axis square
	s = 300; % when we start for the gain analysis
	z = length(data(td).ORN); % where we end
	step_size = 2*act/10;
	if step_size < 0.03
		step_size= 0.03;
	else
		step_size = 0.03*floor(step_size/0.03);
	end
	history_lengths = [0:step_size:2*act];
	example_history_length = history_lengths(3);

	clear x
	x.response = data(td).ORN(s:z);
	x.prediction = data(td).LinearFit(s:z);
	x.stimulus = data(td).PID(s:z);
	x.time = data(td).time(s:z);
	x.filter_length = 201;


	if redo_bootstrap
		data(td).LinearFit_p = GainAnalysis3(x,history_lengths,example_history_length,ph);
	else
		GainAnalysis3(x,history_lengths,example_history_length,ph,NaN*history_lengths);
	end
	clear x


	snapnow;
	delete(gcf);


 clear td

save('/local-data/DA-paper/mahmut_data.mat','data','-append')



figure('outerposition',[0 0 1000 600],'PaperUnits','points','PaperSize',[1300 600]); hold on
ORN =data(1).ORN;
LinearFit=data(1).LinearFit;
PID = data(1).PID;

multiplot(data(1).time,ORN,PID,LinearFit);
set(gca,'XLim',[26 39])

PrettyFig;
