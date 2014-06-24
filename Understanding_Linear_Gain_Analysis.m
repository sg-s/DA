% Understanding Linear Gain Analysis

% load data
load('/local-data/DA-paper/data.mat')

%% Do ORNs rapidly modulate gain?
% Pseudo-white-noise analysis of ORN responses involves presenting binary flickering pulses of odor to the ORN and recording their response. If ORNs rapidly modulate gain on the timescale of response, then responses to pulses of odor in the sequence where the stimulus is locally low will be different from responses to pulses of odor in the sequence where the stimulus is locally high. 

%%
% The data looks like this. The following figure shows the valve state, the odor concentration, and the neuron response. The neuron is ab3A, and the odor presented is 1-octen-3-ol diluted to 3x $10^{-3}$ in Paraffin Oil. The correlation time in the valve position is 30ms. 

%%
% This data file is used for the following analysis:
td = 7;

disp(data(td).original_name)

% detrend PID with a quadratic term
ptrend = fit(data(td).time(:),data(td).PID(:),'Poly2'); 
data(td).PID = data(td).PID(:) - (ptrend(data(td).time) - mean(ptrend(data(td).time)));


% build a simple linear model
[K,~,filtertime] = FindBestFilter(data(td).PID(500:end),data(td).ORN(500:end),[],'filter_length=201;');
data(td).K = K;
data(td).filtertime = filtertime*mean(diff(data(td).time));
data(td).LinearFit = mean(data(td).ORN) + convolve(data(td).time,data(td).PID,data(td).K,data(td).filtertime);


fh=figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(2,1,1), hold on
plot(data(td).time,data(td).PID,'k');
set(gca,'XLim',[min(data(td).time) max(data(td).time)])
ylabel('PID Voltage (V)')
hold off

subplot(2,1,2), hold on
plot(data(td).time,data(td).ORN,'k');
plot(data(td).time,data(td).LinearFit,'r');
set(gca,'XLim',[min(data(td).time) max(data(td).time)])
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
legend ORN LinearFit
PrettyFig;
hold off

snapnow;
delete(fh);

fh=figure('outerposition',[0 0 400 400],'PaperUnits','points','PaperSize',[1000 400]); hold on
plot(data(td).filtertime,data(td).K,'k','LineWidth',2)
title('Linear Filter for this data')
xlabel('FitlerLag (s)')
set(gca,'XLim',[min(data(td).filtertime) max(data(td).filtertime)])
PrettyFig;
hold off

snapnow;
delete(fh);



%%
% The distributions of the input to neuron and the neuron responses are shown below on the left. The autocorrelaiton functions of the input and output are shown on the right. If the ORN modulates its gain on a rapid time-scale, it must do so in this case on a time-scale smaller than the autocorrelation time of the stimulus. 

[act] = PlotDataStatistics(data,td);

snapnow;
delete(gcf);


%%
% Does the ORN selectively amplify responses to relatively small stimuli and suppress responses to relatively high stimuli?
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ph(3) = subplot(1,2,1); hold on 
axis square
ph(4) = subplot(1,2,2); hold on

s = 300; % when we start for the gain analysis
z = length(data(td).ORN); % where we end
example_history_length = 0.111;
history_lengths = [0:mean(diff(data(td).time)):0.02 0.021:0.03:2*act];

clear x
x.response = data(td).ORN(s:z);
x.prediction = data(td).LinearFit(s:z);
x.stimulus = data(td).PID(s:z);
x.time = data(td).time(s:z);
x.filter_length = 201;

redo_bootstrap = 0;
if redo_bootstrap
	data(td).LinearFit_p = GainAnalysis3(x,history_lengths,example_history_length,ph);
else
	GainAnalysis3(x,history_lengths,example_history_length,ph,data(td).LinearFit_p);
end
clear x
hold off


snapnow;
delete(gcf);

% save
save('/local-data/DA-paper/data.mat','data','-append')

%%
% The gain as a function of history length is pretty weird. Why does it not go to 1 when history length is 0? Why is there a mysterious bump around the shortest non-zero history length? 

%%
% To figure this out, let's first repeat the gain analysis, but instead of using the ORN response, we use the linear prediction itself. 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ph(3) = subplot(1,2,1); hold on 
axis square
ph(4) = subplot(1,2,2); hold on

s = 300; % when we start for the gain analysis
z = length(data(td).ORN); % where we end
example_history_length = 0.09;
history_lengths = [0:0.03:2*act];

clear x
x.response = data(td).LinearFit(s:z) + randn(length(data(td).LinearFit(s:z)),1);
x.prediction = data(td).LinearFit(s:z);
x.stimulus = data(td).PID(s:z);
x.time = data(td).time(s:z);
x.filter_length = 201;

redo_bootstrap = 1;
if redo_bootstrap
	p_linear_ORN = GainAnalysis3(x,history_lengths,example_history_length,ph);
else
	GainAnalysis3(x,history_lengths,example_history_length,ph,p_linear_ORN);
end
clear x
hold off


snapnow;
delete(gcf);

