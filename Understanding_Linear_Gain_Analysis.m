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
snapnow;
delete(gcf);

%% Synthetic Data: DA Model vs. Linear Gain Analysis
% Do we still see the weird behaviour of the slope plots with synthetic data? Here we generate the synthetic neuron data by using a DA model to simulate neuron output, with parameters roughly chosen to somewhat match actual ORN statistics. 

p.n_z = 2;
p.n_y = 2;
p.C = 0.3; p.A = 9500; p.B = 25; p.r0 = -80;
p.tau_y = 4; p.tau_z = 8;

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
f=DA_integrate2(data(td).PID,p);
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

snapnow;
delete(gcf);


%% Synthetic Data: DA Model with exp. gaussian inputs vs. Linear Gain Analysis
% So we still see this weird thing. Maybe it's because of something in the stimulus? Now we use the same model but now feed it with exponentiated Gaussian inputs. 

p.n_z = 2;
p.n_y = 2;
p.C = 0.3; p.A = 9500; p.B = 25; p.r0 = -80;
p.tau_y = 4; p.tau_z = 8;

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ph(3) = subplot(1,2,1); hold on 
axis square
ph(4) = subplot(1,2,2); hold on


example_history_length = 0.111;
history_lengths = [0:mean(diff(data(td).time)):0.02 0.021:0.03:2*act];

clear x
% generate fake data
stim = exp(rand(3000,1));
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

snapnow;
delete(gcf);


%%
% OK, now the slopes go to unity as the history length goes to 0. What about the stimulus causes the slopes not to go to 1? Perhaps filtering the stimulus makes this effect come back. 

p.n_z = 2;
p.n_y = 2;
p.C = 0.3; p.A = 9500; p.B = 25; p.r0 = -80;
p.tau_y = 4; p.tau_z = 8;

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ph(3) = subplot(1,2,1); hold on 
axis square
ph(4) = subplot(1,2,2); hold on


example_history_length = 0.111;
history_lengths = [0:mean(diff(data(td).time)):0.02 0.021:0.03:2*act];

clear x
% generate fake data
stim_old = exp(rand(3000,1));
stim = filter(ones(1,40)/40,1,stim_old);
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
% set(ph(4),'YLim',[0.7 1.4],'XLim',[-0.01 max(history_lengths)])

snapnow;
delete(gcf);

%%
% Why does this happen? Why do the gain slopes not go to zero when the history length is zero? This makes sense when we consider what the history length actually is. For varying history lengths, we choose varying smoothing windows and smooth the stimulus with this. However, when the history length is 0 in this case, this doesn't mean the stimulus is not being smoothed! Technically, it is being smoothed, already, as the stimulus itself has correlations. 