% LocustEAGAnalysis.m
% analysis of data from Zane Aldworth and Mark Stopfer
% 
% created by Srinivas Gorur-Shandilya at 2:34 , 27 August 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-n

pHeader;


%% Locust EAG Analysis
% In this document we analyze some data from Zane Aldworth and Mark Stopfer where EAG measurements are made from locust antenna while stimulating them with binary odour signals. 

%% Example Data
% This is what the data looks like. In the following figure, the PID signals have been low-passed filtered by a 30ms filter, and the EAG signals have been high-pass filtered above a 2 second timescale. There are five trials in this dataset, and the shading indicates the standard error of the mean. 

% load data
load('/local-data/DA-paper/data-for-paper/fig4/locust/example-data.mat')

% clean up, sub-sample to 1ms
PID = PID1; clear PID1
EAG = EAG1; clear EAG1 

PID = PID(:,1:10:end)';
EAG = EAG(:,1:10:end)';
valve = ODR1(:,1:10:end)';
valve(valve<max(max(valve))/2) = 0;
valve(valve>0) = 1;

% set zero
for i = 1:width(PID)
	PID(:,i) = PID(:,i) - mean(PID(1:300,i));
	EAG(:,i) = EAG(:,i) - mean(EAG(1:300,i));
end

% compute the derivative of the EAG
dEAG = EAG;
for i = 1:width(PID)
	dEAG(:,i) = filtfilt(ones(50,1),50,[0; diff(dEAG(:,i))]);
end

% filter
PID = bandPass(PID,Inf,30);
EAG = bandPass(EAG,2e3,Inf);

t = 1e-3*(1:length(PID));

figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
subplot(2,1,1), hold on
errorShade(t,mean(PID,2),sem(PID'),'Color',[0 0 0]);
ylabel('Stimulus (mV)')

subplot(2,1,2), hold on
errorShade(t,mean(EAG,2),sem(EAG'),'Color',[0 0 0]);
ylabel('\DeltaEAG (mV)')
xlabel('Time (s)')
prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%% 
% We now back out the stimulus to EAG filter using standard methods. The following figure shows the filter extracted from all 5 trials. Shading is the standard error of the mean. The filter is integrating, and is wider than stimulus to ORN filters in the fruit fly. 

offset = 500;

[K, EAG_prediction, gain] = extractFilters(PID,EAG,'filter_length',2e3,'filter_offset',offset);
t = 1e-3*(1:length(K)) - 1e-3*offset;

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
errorShade(t,mean(K,2),sem(K'),'Color',[0 0 0]);
xlabel('Filter Lag (s)')
ylabel('Filter')
prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% The filter has a lot of structure before 0, which is weird. 

%%
% We now use these filters to project the stimulus and make linear predictions of the response. In the following figure, we compare the linear prediction (red) to the actual response in each trial (black). The top row shows the two time traces overlaid, and the bottom row shows the two traces plotted against each other. We can see that the linear prediction does somewhat OK at predicting the EAG response. 

time = 1e-3*(1:length(PID));
figure('outerposition',[0 0 1400 700],'PaperUnits','points','PaperSize',[1400 800]); hold on
for i = 1:5
	subplot(2,5,i), hold on
	plot(time,EAG(:,i),'k')
	plot(time,EAG_prediction(:,i),'r')
	set(gca,'XLim',[5 15])
	title(['Trial # ' oval(i)])
	if i == 1
		ylabel('EAG/prediction')
		xlabel('Time (s)')
	end

	subplot(2,5,i+5), hold on
	plot(EAG_prediction(1e3:10:end,i),EAG(1e3:10:end,i),'k.')
	if i == 1
		ylabel('EAG')
		xlabel('prediction')
	end
	title(['r^2 = ' oval(rsquare(EAG_prediction(1e3:10:end,i),EAG(1e3:10:end,i)))])
end
prettyFig('fs',20);

if being_published
	snapnow
	delete(gcf)
end

%% Timescale of gain control
% We now measure the gain in 100 ms windows over the entire time trace, and see how that correlates with the mean stimulus in some history length. Estimates of gain are computed by fitting a straight line to a plot of response vs. linear prediction over a 100ms window. Gain estimates are accepted only if linear fits has a r-squared value > 0.8.

%%
% In the following figure, I plot the instantaneous gain vs. the mean stimulus over some window size. The window size is indicated in the title of each panel, and the Spearman correlation between the inst. gain and mean stimulus in that window is shown in the legend. 

[inst_gain,e] = findInstGain(mean(EAG,2),mean(EAG_prediction,2),100,10);

history_lengths = [100 400 1600 3200];

inst_gain(e < .8 | inst_gain < 0) = NaN;
inst_gain(1:5e3) = NaN;

figure('outerposition',[0 0 800 801],'PaperUnits','points','PaperSize',[800 801]); hold on
for i = 1:length(history_lengths)
	autoPlot(length(history_lengths),i); hold on
	mean_stim = vectorise(computeSmoothedStimulus(mean(PID,2),history_lengths(i)));
	l = plot(mean_stim,inst_gain,'k.');
	set(gca,'XScale','log','YScale','log')
	x = mean_stim(1:10:end);
	y = inst_gain(1:10:end);
	rm_this = isnan(x) | isnan(y);
	x(rm_this) = []; y(rm_this) = [];
	rho = spear(x,y);
	legend(['\rho= ' oval(rho)])
	title(['History length  = ' oval(history_lengths(i)) ' ms'])
	xlabel('Mean stimulus in window')
	ylabel('Inst. gain ')
end

prettyFig();

if being_published
	snapnow
	delete(gcf)
end


%%
% Now we want to determine the window size when the instantaneous gain is maximally correlated to the mean stimulus. To do this, we repeat this analysis over many window sizes, and plot the Spearman correlation vs. the window size. We also compare the gain to the mean stimulus for the window size where there is best correlation, and see what the relationship between the two is. 

history_lengths = round(logspace(1,log10(3e4),30)); % all the history lengths we look at, in ms
history_lengths(28) = [];
rho = NaN*history_lengths;

for i = 1:length(history_lengths)
	mean_stim = vectorise(computeSmoothedStimulus(mean(PID,2),history_lengths(i)));
	x = mean_stim(1:10:end);
	y = inst_gain(1:10:end);
	rm_this = isnan(x) | isnan(y);
	x(rm_this) = []; y(rm_this) = [];
	rho(i) = spear(x,y);
end

figure('outerposition',[0 0 1001 500],'PaperUnits','points','PaperSize',[1001 500]); hold on
subplot(1,2,1); hold on
plot(history_lengths,rho,'k+')
xlabel('History lengths (ms)')
set(gca,'XScale','log','YScale','linear')
ylabel('\rho (Inst. gain, mean stimulus)')

subplot(1,2,2); hold on
[~,idx] = min(rho);
mean_stim = vectorise(computeSmoothedStimulus(mean(PID,2),history_lengths(idx)));
plot(mean_stim,inst_gain,'k.')
set(gca,'XScale','log','YScale','log')
x = mean_stim(:); y = inst_gain(:);
rm_this = isnan(mean_stim) | isnan(inst_gain);
ff = fit(x(~rm_this),y(~rm_this),'power1','Upper',[Inf -1],'Lower',[-Inf -1]);
l = plot(sort(x(~rm_this)),ff(sort(x(~rm_this))),'r');
legend(l,'Weber-Fechner Law')
xlabel(['Mean Stimulus in ' char(10) oval(history_lengths(idx)) ' ms'])
ylabel('EAG gain (a.u.)')

prettyFig('FixLogX',true);

if being_published
	snapnow
	delete(gcf)
end

%% EAG derivative 
%  Now, we repeat this analysis with the derivative of the EAG. The advantage of using the derivative is we don't have to filter the EAG to remove slow fluctuations (here, we're still smoothing over 50 ms because the derivative is noisy). 

[K2, dEAG_pred, K1_gain] = extractFilters(PID,dEAG,'filter_length',2e3,'filter_offset',offset);

[inst_gain,e] = findInstGain(mean(dEAG,2),mean(dEAG_pred,2),100,10);

inst_gain(e < .8 | inst_gain < 0) = NaN;
inst_gain(1:5e3) = NaN;

history_lengths = round(logspace(1,log10(3e4),30)); % all the history lengths we look at, in ms
history_lengths(28) = [];
rho = NaN*history_lengths;

for i = 1:length(history_lengths)
	mean_stim = vectorise(computeSmoothedStimulus(mean(PID,2),history_lengths(i)));
	x = mean_stim(1:10:end);
	y = inst_gain(1:10:end);
	rm_this = isnan(x) | isnan(y);
	x(rm_this) = []; y(rm_this) = [];
	rho(i) = spear(x,y);
end


figure('outerposition',[0 0 1502 500],'PaperUnits','points','PaperSize',[1502 500]); hold on
subplot(1,3,1); hold on
filtertime = (1:length(K2)) - offset;
plot(filtertime,mean(K2,2),'k')
xlabel('Filter lag (ms)')
ylabel('Filter to EAG derivative')

subplot(1,3,2), hold on
plot(history_lengths,rho,'k+')
xlabel('History lengths (ms)')
set(gca,'XScale','log','YScale','linear')
ylabel('\rho (Inst. gain, mean stimulus)')

subplot(1,3,3), hold on
[~,idx] = min(rho);
mean_stim = vectorise(computeSmoothedStimulus(mean(PID,2),history_lengths(idx)));
plot(mean_stim,inst_gain,'k.')
set(gca,'XScale','log','YScale','log')
x = mean_stim(:); y = inst_gain(:);
rm_this = isnan(mean_stim) | isnan(inst_gain);
ff = fit(x(~rm_this),y(~rm_this),'power1','Upper',[Inf -1],'Lower',[-Inf -1]);
l = plot(sort(x(~rm_this)),ff(sort(x(~rm_this))),'r');
legend(l,'Weber-Fechner Law')
xlabel(['Mean Stimulus in ' char(10) oval(history_lengths(idx)) ' ms'])
ylabel('EAG gain (a.u.)')

prettyFig('FixLogX',true);

if being_published
	snapnow
	delete(gcf)
end

%% Whiff-based analysis
% Now, we perform a whiff-based analysis. Instead of finding the "instantenous" gain, which is hard to measure, and not very instantaneous, we define gains over whiffs of odorants, identified as rising phases in the stimulus. 


history_lengths = round(logspace(1,log10(3e4),30)); % all the history lengths we look at, in ms
history_lengths(28) = [];

pred = nanmean(EAG_prediction,2); 
resp = nanmean(EAG,2);  
stim = nanmean(PID,2); 

% find when the valve opens
[ons,offs] = findWhiffs(stim);
ons = ons+finddelay(stim,resp);
offs = offs+finddelay(stim,resp);

rm_this = ons>length(stim);
ons(rm_this) = [];
offs(rm_this) = [];

% plot the gain in each of these windows
[gain,gain_err] = findGainInWindows(ons,offs,pred,resp);

rm_this = gain < 0 | gain_err < .5;
gain(rm_this) = [];
ons(rm_this) = [];
offs(rm_this) = [];
gain_err(rm_this) = [];

% find the mean stimulus in the preceding X ms in these windows

% also find rho for various values of the history length and plot it
rho = findRhoForHistoryLengths(gain,stim,ons,offs,history_lengths);

figure('outerposition',[0 0 501 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(history_lengths,rho,'k+')
xlabel('history length (ms)')
title('Gain computed during "whiffs"')
ylabel('\rho')
set(gca,'XScale','log')

prettyFig('FixLogX',true);

if being_published
	snapnow
	delete(gcf)
end


%% Version Info
% The file that generated this document is called:
disp(mfilename)

%%
% and its md5 hash is:
Opt.Input = 'file';
disp(dataHash(strcat(mfilename,'.m'),Opt))

%%
% This file should be in this commit:
[status,m]=unix('git rev-parse HEAD');
if ~status
	disp(m)
end

%%
% This file has the following external dependencies:
showDependencyHash(mfilename);

t = toc;

%% 
% This document was built in: 
disp(strcat(oval(t,3),' seconds.'))

% tag the file as being published 

if being_published
	unix(['tag -a published ',which(mfilename)]);
	unix(['tag -r publish-failed ',which(mfilename)]);
end
