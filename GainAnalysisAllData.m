% GainAnalysisAllData.m
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


% some parameters
font_size = 20;
marker_size = 10;
marker_size2 = 20;

% which data set to use?
td = 3;

redo_bootstrap = 1;

% load data
load('/local-data/DA-paper/data.mat')

%                     ########     ###    ########    ###    
%                     ##     ##   ## ##      ##      ## ##   
%                     ##     ##  ##   ##     ##     ##   ##  
%                     ##     ## ##     ##    ##    ##     ## 
%                     ##     ## #########    ##    ######### 
%                     ##     ## ##     ##    ##    ##     ## 
%                     ########  ##     ##    ##    ##     ## 


%% Data Overview
% The data looks like this. The following figure shows the valve state, the odor concentration, and the neuron response. The stimulus measurement is de-trended with a quadratic term to remove what is supposed to be a drift in the PID. 

%%
% This data file is used for the following analysis:


disp(data(td).original_name)

% detrend PID with a quadratic term
ptrend = fit(data(td).time(:),data(td).PID(:),'Poly2'); 
data(td).PID = data(td).PID(:) - (ptrend(data(td).time) - mean(ptrend(data(td).time)));


fh=figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(2,1,1), hold on
plot(data(td).time,data(td).PID,'k');
set(gca,'XLim',[min(data(td).time) max(data(td).time)])
ylabel('PID Voltage (V)')
hold off

subplot(2,1,2), hold on
plot(data(td).time,data(td).ORN,'k');
set(gca,'XLim',[min(data(td).time) max(data(td).time)])
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
PrettyFig;
hold off

snapnow;
delete(fh);


%      ######  ########    ###    ######## ####  ######  ######## ####  ######   ######  
%     ##    ##    ##      ## ##      ##     ##  ##    ##    ##     ##  ##    ## ##    ## 
%     ##          ##     ##   ##     ##     ##  ##          ##     ##  ##       ##       
%      ######     ##    ##     ##    ##     ##   ######     ##     ##  ##        ######  
%           ##    ##    #########    ##     ##        ##    ##     ##  ##             ## 
%     ##    ##    ##    ##     ##    ##     ##  ##    ##    ##     ##  ##    ## ##    ## 
%      ######     ##    ##     ##    ##    ####  ######     ##    ####  ######   ######  

%% Data Statistics
% In the following section, we plot some statistical features of the stimulus. In the figure below, the panel on the left shows the distribution (normalised) of the values of the stimulus and the response. The panel on the right shows the autocorrelation function of the stimulus and the response. From the autocorrelation function, we can also estimate the time at which the autocorrelation of the response is zero, which we will use to in future analyses.

clear ph
figure('outerposition',[0 0 1100 500],'PaperUnits','points','PaperSize',[1100 500]); hold on
ph(1)=subplot(1,2,1); hold on
ph(2)=subplot(1,2,2); hold on
act = PlotDataStatistics(data,td,ph);

set(ph(2),'XScale','linear')

PrettyFig;

snapnow;
delete(gcf);


%      ##       #### ##    ## ########    ###    ########           ######## #### ######## 
%      ##        ##  ###   ## ##         ## ##   ##     ##          ##        ##     ##    
%      ##        ##  ####  ## ##        ##   ##  ##     ##          ##        ##     ##    
%      ##        ##  ## ## ## ######   ##     ## ########           ######    ##     ##    
%      ##        ##  ##  #### ##       ######### ##   ##            ##        ##     ##    
%      ##        ##  ##   ### ##       ##     ## ##    ##           ##        ##     ##    
%      ######## #### ##    ## ######## ##     ## ##     ##          ##       ####    ##    


%% Linear Fit to Data
% We can back out a filter from this data, and it looks like this:


% build a simple linear model
[K,~,filtertime] = FindBestFilter(data(td).PID(500:end),data(td).ORN(500:end),[],'filter_length=201;');
data(td).K = K;
data(td).filtertime = filtertime*mean(diff(data(td).time));
data(td).LinearFit = convolve(data(td).time,data(td).PID,data(td).K,data(td).filtertime);


fh=figure('outerposition',[0 0 400 400],'PaperUnits','points','PaperSize',[1000 400]); hold on
plot(data(td).filtertime,data(td).K,'k','LineWidth',2)
title('Linear Filter for this data')
xlabel('FitlerLag (s)')
ylabel('Filter Amplitude (Hz)')
set(gca,'XLim',[min(data(td).filtertime) max(data(td).filtertime)])

PrettyFig;

snapnow;
delete(fh);

%%
% What does the prediction of this filter look like? The following figure shows the data (black) and the linear prediction (red) superimposed on it. 

figure('outerposition',[0 0 1100 400],'PaperUnits','points','PaperSize',[1100 400]); hold on
plot(data(td).time,data(td).ORN,'k');
plot(data(td).time,data(td).LinearFit,'r');
set(gca,'XLim',[min(data(td).time) max(data(td).time)])
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')

PrettyFig;

snapnow;
delete(fh);

%%
% Once again, we need to correct the linear prediction by adding back the mean of the response. We also pass it through a rectifier, because negative values have no meaning.  

data(td).LinearFit = data(td).LinearFit + mean(data(td).ORN);
data(td).LinearFit(data(td).LinearFit<0)=0;

figure('outerposition',[0 0 1100 400],'PaperUnits','points','PaperSize',[1100 400]); hold on
plot(data(td).time,data(td).ORN,'k');
plot(data(td).time,data(td).LinearFit,'r');
set(gca,'XLim',[min(data(td).time) max(data(td).time)])
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')

PrettyFig;

snapnow;
delete(fh);


%%
% The r-square of the linear prediction is: 

disp(rsquare(data(td).LinearFit,data(td).ORN))

%%
% And the raw Euclidian Distance between the prediction and the data is: 

disp(Cost2(data(td).LinearFit(205:end-33),data(td).ORN(205:end-33)))


%         #######  ##     ## ######## ########  ##     ## ########    ##    ## ##       
%        ##     ## ##     ##    ##    ##     ## ##     ##    ##       ###   ## ##       
%        ##     ## ##     ##    ##    ##     ## ##     ##    ##       ####  ## ##       
%        ##     ## ##     ##    ##    ########  ##     ##    ##       ## ## ## ##       
%        ##     ## ##     ##    ##    ##        ##     ##    ##       ##  #### ##       
%        ##     ## ##     ##    ##    ##        ##     ##    ##       ##   ### ##       
%         #######   #######     ##    ##         #######     ##       ##    ## ######## 

%% Output Non-linearity
% Does adding an output non-linearity post-hoc improve the linear fit? In the following section, we fit a two-parameter Hill function to the linear prediction and the data, after we have generated the linear prediction. (In other words, the Hill function is fit _after_ the best linear prediction is calculated.) The following figure shows the shape of the best-fit non-linearity:

xdata = data(td).LinearFit;
ydata = data(td).ORN;

% crop it to lose NaNs
ydata(isnan(xdata)) = [];
xdata(isnan(xdata)) = [];

xdata = xdata(:);
ydata = ydata(:);

% rectify
xdata(xdata<0)=0;

fo=optimset('MaxFunEvals',1000,'Display','none');
x = lsqcurvefit(@hill,[50 2 2],xdata,ydata,[],[],fo);

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(xdata,hill(x,xdata),'k')
xlabel('Linear Prediction (Hz)')
ylabel('Nonlinearity Output (Hz)')

PrettyFig;

snapnow;
delete(gcf);

%%
% How does the output nonlinearity change the prediction? In the following figure, the data is shown in black, the linear prediction is shown in red, and the LN prediction is shown in green. 

% save this for later
xdata = data(td).LinearFit;
xdata(xdata<0)=0;
LNpred = hill(x,xdata);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(data(td).time,data(td).ORN,'k')
plot(data(td).time,data(td).LinearFit,'r')
plot(data(td).time,hill(x,LNpred),'g')
set(gca,'XLim',[18 22])
xlabel('Time (s)')
ylabel('Firing rate (Hz)')

PrettyFig;

snapnow;
delete(gcf);



%%
%The r-square of the LN prediction is: 

disp(rsquare(LNpred,data(td).ORN))

%%
% And the raw Euclidian Distance between the prediction and the data is: 

disp(Cost2(LNpred(205:end-33),data(td).ORN(205:end-33)))


%         ######      ###    #### ##    ##       ###    ##    ##    ###    ##       
%        ##    ##    ## ##    ##  ###   ##      ## ##   ###   ##   ## ##   ##       
%        ##         ##   ##   ##  ####  ##     ##   ##  ####  ##  ##   ##  ##       
%        ##   #### ##     ##  ##  ## ## ##    ##     ## ## ## ## ##     ## ##       
%        ##    ##  #########  ##  ##  ####    ######### ##  #### ######### ##       
%        ##    ##  ##     ##  ##  ##   ###    ##     ## ##   ### ##     ## ##       
%         ######   ##     ## #### ##    ##    ##     ## ##    ## ##     ## ######## 

%% Linear Model Gain Analysis
% In this section, we perform a gain analysis on the linear prediction first for an arbitrarily chosen history length of 120ms (panel on the left) and then for various history lengths (panel on the right). The history lengths where the slopes are significantly different (p<0.05) are indicated by dots. Significance is determined by bootstrapping the data 100 times. History lengths up to twice the autocorrelation length of the stimulus are investigated. 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ph(3) = subplot(1,2,1); hold on 
axis square
ph(4) = subplot(1,2,2); hold on

s = 300; % when we start for the gain analysis
z = length(data(td).ORN) - 33; % where we end
example_history_length = 0.12;
history_lengths = [0:0.06:2*act];

clear x
x.response = data(td).ORN(s:z);
x.prediction = data(td).LinearFit(s:z);
x.stimulus = data(td).PID(s:z);
x.time = data(td).time(s:z);
x.filter_length = 201;


if redo_bootstrap
	ptemp = GainAnalysis3(x,history_lengths,example_history_length,ph);
else
	GainAnalysis3(x,history_lengths,example_history_length,ph,ptemp);
end
clear x


snapnow;
delete(gcf);

%       ########  ##     ## ##        ######  ########       ###    ##     ## ########  
%       ##     ## ##     ## ##       ##    ## ##            ## ##   ###   ### ##     ## 
%       ##     ## ##     ## ##       ##       ##           ##   ##  #### #### ##     ## 
%       ########  ##     ## ##        ######  ######      ##     ## ## ### ## ########  
%       ##        ##     ## ##             ## ##          ######### ##     ## ##        
%       ##        ##     ## ##       ##    ## ##          ##     ## ##     ## ##        
%       ##         #######  ########  ######  ########    ##     ## ##     ## ##        

%% Gain Analysis using Pulse Amplitudes
% This doesn't look very good, because this method of gain analysis doesn't work well when the response of the neuron (and the prediction) is uniformly zero. In this section, we check for fast gain control using a different method. In the following figure, the amplitude of the peak of each pulse (in response to odor) is plotted vs. the length of time preceding the pulse during which the valve was off (implying no stimulus was being delivered to the neuron.)

%%
% This correlation is not due to the stimulus being different, as can be seen from the panel on the right. 

[ons,offs] = ComputeOnsOffs(data(td).Valve);
PulsePeaks = NaN(1,length(ons));
PIDPeaks = NaN(1,length(ons));


a = ons; b = offs;
a(1) = []; b(end) = [];
ValveOffDurations = (a-b)*mean(diff(data(td).time));

for i = 2:length(ons)
	PulsePeaks(i)=max(data(td).ORN(ons(i):offs(i)));
	PIDPeaks(i)=max(data(td).PID(ons(i):offs(i)));
end
PulsePeaks(1) = [];
PIDPeaks(1) = [];

% fit a line
[ffit, gof] = fit(ValveOffDurations(:),PulsePeaks(:),'Poly1');
x = sort(unique(ValveOffDurations));
y = ffit(x);


figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
plot(ValveOffDurations,PulsePeaks,'k.','MarkerSize',24)
plot(x,y,'r')
ylabel('Peak response (Hz)')
xlabel('Preceding Blank Duration (s)')


subplot(1,2,2), hold on
plot(ValveOffDurations,PIDPeaks,'k.','MarkerSize',24)
ylabel('Peak stimulus of pulse (a.u.)')
xlabel('Preceding Blank Duration (s)')

PrettyFig;

snapnow;
delete(gcf);


%%
% and the r-square of the fit to response height vs. blank duration is:
disp(gof.rsquare)

%%
% In contrast, the r-square of the fit to stimulus peak vs. blank duration is:
disp(rsquare(ValveOffDurations,PIDPeaks))

return

% ##       ##    ##     ######      ###    #### ##    ##       ###    ##    ##    ###    ##       
% ##       ###   ##    ##    ##    ## ##    ##  ###   ##      ## ##   ###   ##   ## ##   ##       
% ##       ####  ##    ##         ##   ##   ##  ####  ##     ##   ##  ####  ##  ##   ##  ##       
% ##       ## ## ##    ##   #### ##     ##  ##  ## ## ##    ##     ## ## ## ## ##     ## ##       
% ##       ##  ####    ##    ##  #########  ##  ##  ####    ######### ##  #### ######### ##       
% ##       ##   ###    ##    ##  ##     ##  ##  ##   ###    ##     ## ##   ### ##     ## ##       
% ######## ##    ##     ######   ##     ## #### ##    ##    ##     ## ##    ## ##     ## ######## 


%% LN Model Gain Analysis
% In the previous section, we saw that the simple linear model fails to account for this fast adaptation of the ORNs. Specifically, the gain of the neuron w.r.t to the model is significantly different for times when the stimulus is high and the when the stimulus is low.

%%
% In this section, we want to know if the LN model (adding a non-linear function post-hoc) corrects for this misprediction of gain. Here, we repeat the gain analysis as in the previous section, but this time, using the LN prediction instead of the linear prediction. 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ph(3) = subplot(1,2,1); hold on 
axis square
ph(4) = subplot(1,2,2); hold on

s = 300; % when we start for the gain analysis
z = length(data(td).ORN) - 33; % where we end
example_history_length = 0.12;
history_lengths = [0:0.06:2*act];

clear x

x.response = data(td).ORN(s:z);
x.prediction = LNpred(s:z);
x.stimulus = data(td).PID(s:z);
x.time = data(td).time(s:z);
x.filter_length = 201;


if redo_bootstrap
	ptemp2 = GainAnalysis3(x,history_lengths,example_history_length,ph);
else
	GainAnalysis3(x,history_lengths,example_history_length,ph,ptemp2);
end
clear x


snapnow;
delete(gcf);



%% Understanding features of Linear Model Gain Analysis
% There are some features of the gain analysis that are interesting, especially when compared to the gain analysis of the NL model: 
% 
% * The slopes of the linear model gain analysis do not go to 1 even when the history lengths go to 0
% * The LN model does not behave this way: it behaves as expected, and the slopes go to zero as the history length go to 0
% * The green and the red curves are not mirror images of each other. Why not? What makes them different? 

%% 
% In the following figure, we use a DA model to generate synthetic responses to the actual stimulus and then fit a linear filter to this data, and then perform the gain analysis as above on this synthetic dataset. 

%% 
% The parameters of the DA model are chosen to be a best fit to the actual experimental data, so the DA model chosen best represents the actual neuron's response. 

if ~exist('p')
	d.stimulus = data(td).PID;
	d.response = data(td).ORN;
	[p, f] = FitDAModelToData(d,[25576 85 0.03 1.5 2 13 2 -151],[900 50 0 1e-1 2 1e-2 2 -200],[98000 5000 1 20 2 40 2 1],200);
end


figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ph(3) = subplot(1,2,1); hold on 
axis square
ph(4) = subplot(1,2,2); hold on

s = 300; % when we start for the gain analysis
z = length(data(td).ORN); % where we end
example_history_length = 0.12;
history_lengths = [0:0.06:2*act];

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
	GainAnalysis3(x,history_lengths,example_history_length,ph,p_fake_ORN);
end
clear x
hold off
set(ph(4),'YLim',[0.7 1.4],'XLim',[-0.01 max(history_lengths)])
title(ph(4),'DA Model response to odor stimulus')
ylabel(ph(3),'DA Model Response')

snapnow;
delete(gcf);

%%
% The slopes still don't go to 1, and now there is a weird deflection at low history lengths. What if we shuffle the stimulus? 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ph(3) = subplot(1,2,1); hold on 
axis square
ph(4) = subplot(1,2,2); hold on

s = 300; % when we start for the gain analysis
z = length(data(td).ORN); % where we end
example_history_length = 0.12;
history_lengths = [0:0.06:2*act];

clear x
% generate fake data
f=DA_integrate2(d.stimulus,p);
[K,~,filtertime] = FindBestFilter(data(td).PID,f);
fp=mean(f(s:z))+convolve(data(td).time,data(td).PID,K,filtertime);



x.response = f(s:z);
x.prediction = fp(s:z);
x.stimulus = data(td).PID(s:z);
x.stimulus=x.stimulus(randperm(length(x.stimulus)));
x.time = data(td).time(s:z);
x.filter_length = length(K);

if redo_bootstrap
	p_fake_ORN2 = GainAnalysis3(x,history_lengths,example_history_length,ph);
else
	GainAnalysis3(x,history_lengths,example_history_length,ph,p_fake_ORN2);
end
clear x
hold off
set(ph(4),'YLim',[0.7 1.4],'XLim',[-0.01 max(history_lengths)])
title(ph(4),'DA Model response to shuffled stimulus')
ylabel(ph(3),'DA Model Response')

snapnow;
delete(gcf);

%%
% Now, the slopes go to 1, but then, the slopes are always 1 no matter what the history length. We want to create a condition where the slopes. What if we shuffle, then filter the odor stimulus? Perhaps correlation times in the stimulus affect this slope behaviour? The following figure shows the effect of shuffling and filtering the stimulus:

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
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


%%
% We now repeat the gain analysis on this new synthetic dataset. 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ph(3) = subplot(1,2,1); hold on 
axis square
ph(4) = subplot(1,2,2); hold on

s = 300; % when we start for the gain analysis
z = length(data(td).ORN); % where we end
example_history_length = 0.12;
history_lengths = [0:0.06:2*act];

clear x
% generate fake data
f=DA_integrate2(d.stimulus,p);
[K,~,filtertime] = FindBestFilter(data(td).PID,f);
fp=mean(f(s:z))+convolve(data(td).time,data(td).PID,K,filtertime);




x.response = f(s:z);
x.prediction = fp(s:z);
x.stimulus = stim(s:z);
x.time = data(td).time(s:z);
x.filter_length = length(K);

if redo_bootstrap
	p_fake_ORN3 = GainAnalysis3(x,history_lengths,example_history_length,ph);
else
	GainAnalysis3(x,history_lengths,example_history_length,ph,p_fake_ORN3);
end
clear x
hold off
set(ph(4),'YLim',[0.7 1.4],'XLim',[-0.01 max(history_lengths)])
title(ph(4),'Response to shuffled filtered stimulus')
ylabel(ph(3),'DA Model Response')

snapnow;
delete(gcf);
