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
redo_bootstrap = 1;

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

if redo_bootstrap
	data(td).LinearFit_p = GainAnalysis3(x,history_lengths,example_history_length,ph);
else
	GainAnalysis3(x,history_lengths,example_history_length,ph,Inf*history_lengths);
end
clear x
hold off
set(ph(4),'YLim',[0.7 1.4],'XLim',[-0.01 max(history_lengths)])

title(ph(4),'ORN response to odor stimulus')

snapnow;
delete(gcf);

% save
save('/local-data/DA-paper/data.mat','data','-append')

%%
% The gain as a function of history length is pretty weird. Why does it not go to 1 when history length is 0? Why is there a mysterious downward deflection close to 0?

%% Sanity check: Linear synthetic data vs. linear gain analysis
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
x.response = data(td).LinearFit(s:z);
x.prediction = data(td).LinearFit(s:z);
x.stimulus = data(td).PID(s:z);
x.time = data(td).time(s:z);
x.filter_length = 201;

if redo_bootstrap
	p_linear_ORN = GainAnalysis3(x,history_lengths,example_history_length,ph);
else
	GainAnalysis3(x,history_lengths,example_history_length,ph,Inf*history_lengths);
end
clear x
hold off
set(ph(4),'YLim',[0.7 1.4],'XLim',[-0.01 max(history_lengths)])
title(ph(4),'Linear Model Prediction to odor stimulus')
ylabel(ph(3),'Linear Model Prediction')

snapnow;
delete(gcf);

%% Synthetic Data: DA Model vs. Linear Gain Analysis
% Do we still see the weird behaviour of the slope plots with synthetic data? Here we generate the synthetic neuron data by using a DA model to simulate neuron output, with parameters derived from a fit of the DA model to the data. The parameters of the DA model are such that the output fits the data quite well, the r-square to the data is:

disp(rsquare(Rguess(300:end),data(td).ORN(300:end)))

if ~exist('p')
	d.stimulus = data(td).PID;
	d.response = data(td).ORN;
	[p, f] = FitDAModelToData(d,[2.6422e+04 90 0.03 1.5 2 13 2 -151],[900 50 0 1e-1 1 1e-2 1 -200],[98000 5000 1 20 10 40 10 1],200);
end


figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ph(3) = subplot(1,2,1); hold on 
axis square
ph(4) = subplot(1,2,2); hold on

s = 300; % when we start for the gain analysis
z = length(data(td).ORN); % where we end
example_history_length = 0.111;
history_lengths = [0:mean(diff(data(td).time)):0.02 0.021:0.03:2*act];

clear x
% generate fake data
f=DA_integrate2(d.stimulus,p);
[K,~,filtertime] = FindBestFilter(data(td).PID,f);
fp=mean(f(s:z))+convolve(data(td).time,data(td).PID,K,filtertime);


x.response = f(s:z);
x.prediction = fp(s:z);
x.stimulus = data(td).PID(s:z);
x.time = data(td).time(s:z);
x.filter_length = length(K);

if redo_bootstrap
	p_fake_ORN = GainAnalysis3(x,history_lengths,example_history_length,ph);
else
	GainAnalysis3(x,history_lengths,example_history_length,ph,Inf*history_lengths);
end
clear x
hold off
set(ph(4),'YLim',[0.7 1.4],'XLim',[-0.01 max(history_lengths)])
title(ph(4),'DA Model response to odor stimulus')
ylabel(ph(3),'DA Model Response')

snapnow;
delete(gcf);

%%
% So we see that the DA model also captures, qualitatively, the two interesting things in the real data: 1) the fact that the lines do not go to 1 as the history lengths go to zero and 2) the weird downward deflection close to zero. 


%% Synthetic Data: DA Model with exp. gaussian inputs vs. Linear Gain Analysis
% Maybe it's because of something in the stimulus? Now we use the same DA model but now feed it with exponentiated Gaussian inputs. 


figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ph(3) = subplot(1,2,1); hold on 
axis square
ph(4) = subplot(1,2,2); hold on


example_history_length = 0.111;
history_lengths = [0:mean(diff(data(td).time)):0.02 0.021:0.03:2*act];

clear x
% generate fake data
stim = exp(rand(3000,1))/50;
t = mean(diff(data(td).time)):mean(diff(data(td).time)):mean(diff(data(td).time))*length(stim);
f=DA_integrate2(stim,p);
[K,~,filtertime] = FindBestFilter(stim,f,[],'filter_length=200;');

s = 300; % when we start for the gain analysis
z = length(stim); % where we end
fp=mean(f(s:z))+convolve(t,stim,K,filtertime);

x.response = f(s:z);
x.prediction = fp(s:z);
x.stimulus = stim(s:z);
x.time = t(s:z);
x.filter_length = length(K);

if redo_bootstrap
	p_fake_ORN2 = GainAnalysis3(x,history_lengths,example_history_length,ph);
else
	GainAnalysis3(x,history_lengths,example_history_length,ph,Inf*history_lengths);
end
clear x
hold off
set(ph(4),'YLim',[0.7 1.4],'XLim',[-0.01 max(history_lengths)])

title(ph(4),'DA Model response to uncorrelated noise')
ylabel(ph(3),'DA Model Response')

snapnow;
delete(gcf);


%%
% OK, now the slopes go to unity as the history length goes to 0. What about the stimulus causes the slopes not to go to 1? Here we present a different stimulus to the DA model to see how the gain analysis reports the gain at short time scales. The stimulus presented is the actual PID trace, shuffled, and rescaled so that the magnitudes of ORN responses are approximately the same. 



figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ph(3) = subplot(1,2,1); hold on 
axis square
ph(4) = subplot(1,2,2); hold on


example_history_length = 0.111;
history_lengths = [0:mean(diff(data(td).time)):0.02 0.021:0.03:2*act];

clear x
% generate fake data
stim=d.stimulus(randperm(length(d.stimulus)));
stim = stim-min(stim)/2;
stim = stim*2;
t = mean(diff(data(td).time)):mean(diff(data(td).time)):mean(diff(data(td).time))*length(stim);
f=DA_integrate2(stim,p);

f(1:100) = [];
t(1:100) = [];
stim(1:100) = [];

s = 300; % when we start for the gain analysis
z = length(stim); % where we end

[K,~,filtertime] = FindBestFilter(stim(s:z),f(s:z),[],'filter_length=200;');
fp=mean(f)+convolve(t,stim,K,filtertime);



x.response = f(s:z);
x.prediction = fp(s:z);
x.stimulus = stim(s:z);
x.time = t(s:z);
x.filter_length = length(K);

if redo_bootstrap
	p_fake_ORN2 = GainAnalysis3(x,history_lengths,example_history_length,ph);
else
	GainAnalysis3(x,history_lengths,example_history_length,ph,Inf*history_lengths);
end
clear x
hold off
set(ph(4),'YLim',[0.7 1.4],'XLim',[-0.01 max(history_lengths)])
title(ph(4),'DA Model response to shuffed odor')
ylabel(ph(3),'DA Model Response')

snapnow;
delete(gcf);

%%
% OK, the gain analysis still reports that the gains go to 1 as the history length goes to zero. What manipulation can break this, and recapitulate the effect we see with the real data? Perhaps filtering will do the trick. In the following plots, the stimulus being presented to the DA model is the actual PID trace, shuffled, and filtered, and then rescaled so that it matches the real PID closely in both amplitude and temporal structure. 



figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1); hold on 
stim=d.stimulus(randperm(length(d.stimulus)));
stim = filter(ones(100,1)/100,1,stim-mean(stim))+mean(stim);
stim = stim-mean(stim);
stim = stim*10;
stim = stim+mean(d.stimulus);
plot(d.stimulus,'k'), hold on
plot(stim,'r')
set(gca,'XLim',[2000 4000])
xlabel('Time (a.u.)')
ylabel('Stimulus (a.u.)')
legend PID ShuffledFilteredPID
PrettyFig;

subplot(1,2,2); hold on 
[y] = autocorr(d.stimulus,150);
atime = 0:mean(diff(data(td).time)):mean(diff(data(td).time))*150;
plot(atime,y,'k')
[y] = autocorr(stim,150);
atime = 0:mean(diff(data(td).time)):mean(diff(data(td).time))*150;
plot(atime,y,'r')
legend PID ShuffledFilteredPID
set(gca,'XLim',[0 0.5])
xlabel('Lag (s)')
ylabel('Autocorrelation')
PrettyFig;


figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ph(3) = subplot(1,2,1); hold on 
axis square
ph(4) = subplot(1,2,2); hold on


example_history_length = 0.111;
history_lengths = [0:mean(diff(data(td).time)):0.02 0.021:0.03:2*act];

clear x
t = mean(diff(data(td).time)):mean(diff(data(td).time)):mean(diff(data(td).time))*length(stim);
f=DA_integrate2(stim,p);

f(1:100) = [];
t(1:100) = [];
stim(1:100) = [];

s = 300; % when we start for the gain analysis
z = length(stim); % where we end

[K,~,filtertime] = FindBestFilter(stim(s:z),f(s:z),[],'filter_length=200;');
fp=mean(f)+convolve(t,stim,K,filtertime);



x.response = f(s:z);
x.prediction = fp(s:z);
x.stimulus = stim(s:z);
x.time = t(s:z);
x.filter_length = length(K);

if redo_bootstrap
	p_fake_ORN2 = GainAnalysis3(x,history_lengths,example_history_length,ph);
else
	GainAnalysis3(x,history_lengths,example_history_length,ph,Inf*history_lengths);
end
clear x
hold off
set(ph(4),'YLim',[0.7 1.4],'XLim',[-0.01 max(history_lengths)])
title(ph(4),'DA Model response to shuffled, filtered odor stimulus')
ylabel(ph(3),'DA Model Response')

snapnow;
delete(gcf);
