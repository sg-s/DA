% GainAnalysisAllData.m
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

%%
% The purpose of this document is to determine if ORNs show the characteristic features of fast gain control in flickering stimuli experiments with fast odors, and with other receptors. We want to answer the question of how general the phenomenon of fast gain control that we saw with ab3A and 1-octen-3-ol is. 

% some parameters
font_size = 20;
marker_size = 10;
marker_size2 = 20;

% which data set to use?
td = 5;

% define history lengths to perform gain analysis on
history_lengths=[0:0.03:0.3 0.36:0.06:1 1.2:1.2:5];

% define this to either recompute everything, or set to zero if you're ready to publish
% this is commented out so that you have to explicitly specify it so you know what you're doing.
% redo_bootstrap = 1;


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
% How does the output non-linearity change the prediction? In the following figure, the data is shown in black, the linear prediction is shown in red, and the LN prediction is shown in green. 

% save this for later
xdata = data(td).LinearFit;
xdata(xdata<0)=0;
LNpred = hill(x,xdata);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(data(td).time,data(td).ORN,'k')
plot(data(td).time,data(td).LinearFit,'r')
plot(data(td).time,hill(x,LNpred),'g')
set(gca,'XLim',[19 23])
xlabel('Time (s)')
ylabel('Firing rate (Hz)')

PrettyFig;

snapnow;
delete(gcf);



%%
%The r-square of the LN prediction is: 

disp(rsquare(LNpred,data(td).ORN))

%%
% And the raw Euclidean Distance between the prediction and the data is: 

disp(Cost2(LNpred(205:end-33),data(td).ORN(205:end-33)))


%         ######      ###    #### ##    ##       ###    ##    ##    ###    ##       
%        ##    ##    ## ##    ##  ###   ##      ## ##   ###   ##   ## ##   ##       
%        ##         ##   ##   ##  ####  ##     ##   ##  ####  ##  ##   ##  ##       
%        ##   #### ##     ##  ##  ## ## ##    ##     ## ## ## ## ##     ## ##       
%        ##    ##  #########  ##  ##  ####    ######### ##  #### ######### ##       
%        ##    ##  ##     ##  ##  ##   ###    ##     ## ##   ### ##     ## ##       
%         ######   ##     ## #### ##    ##    ##     ## ##    ## ##     ## ######## 

%% Linear Model Gain Analysis
% In this section, we perform a gain analysis on the linear prediction. 

%%
% First figure: The top panel shows the stimulus (black), and the smoothed stimulus smoothed over the history window (used in the example in the second figure). The top10% of the smoothed stimulus are highlighted in red, and the bottom 10% are highlighted in green. For these time points, the data (black, bottom panel) and the the prediction (red, bottom panel) are compared. 

%%
% Second figure: In the left panel, these the prediction and the data at these chosen time points are plotted, along with all the data in gray, together with best fit linear functions. The panel on the right shows the same analysis for all history lengths, with dots indicating points where the slopes between the red and green clouds are significantly different (p<0.05, data bootstrapped a 100 times). 

f1=figure('outerposition',[0 0 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
ph(1) = subplot(2,1,1); hold on 
ph(2) = subplot(2,1,2); hold on

title(ph(1),strrep(data(td).original_name,'_','-'),'FontSize',20);

f2=figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ph(3) = subplot(1,2,1); hold on 
axis square
ph(4) = subplot(1,2,2); hold on



s = 300; % when we start for the gain analysis
z = length(data(td).ORN) - 33; % where we end
example_history_length = 0.96;




clear x
x.response = data(td).ORN(s:z);
x.prediction = LNpred(s:z);
x.stimulus = data(td).PID(s:z);
x.time = data(td).time(s:z);
x.filter_length = 201;


if redo_bootstrap
	ptemp = GainAnalysis3(x,history_lengths,example_history_length,ph);
else
	GainAnalysis3(x,history_lengths,example_history_length,ph,ptemp);
end

snapnow;
delete(f1);

snapnow;
delete(f2);




%%
% There are many things different between A) the example 1-octen-3-ol dataset and B) the current dataset.
% 
% * In A, ORN firing rates never go to 0. In B, they do.
% * In A, ORN firing rates when the stimulus is in lowest 10% is strictly > 0. Not so in B.
% * The stimulus distribution is unimodal in A, with fluctuations around a non-zero mean. In B, the stimulus distribution is sharply bimodal, with the stimulus either being zero or close to the maximum. 
% * In A, ORN firing rates with the stimulus in the highest 10% is close to the maximum. In B, it is far from the maximum, as the neuron has adapted to the onset of the pulse. 
% * In B, the gain of the filter is not exactly 1 as the filter makes negative predictions, which we set to 0 as negative firing rates don't make sense. 


%       ########  ##     ## ##        ######  ########       ###    ##     ## ########  
%       ##     ## ##     ## ##       ##    ## ##            ## ##   ###   ### ##     ## 
%       ##     ## ##     ## ##       ##       ##           ##   ##  #### #### ##     ## 
%       ########  ##     ## ##        ######  ######      ##     ## ## ### ## ########  
%       ##        ##     ## ##             ## ##          ######### ##     ## ##        
%       ##        ##     ## ##       ##    ## ##          ##     ## ##     ## ##        
%       ##         #######  ########  ######  ########    ##     ## ##     ## ##        

%% Gain Analysis using Pulse Amplitudes
% In this section, we check for fast gain control using a different method. In the following figure, the amplitude of the peak of each pulse (in response to odor) is plotted vs. the length of time preceding the pulse during which the valve was off (implying no stimulus was being delivered to the neuron.)

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

%%
% The plots above show that the longer the low stimulus period is preceding the current pulse, the large is the response of the neuron, given a constant pulse height. The duration of the preceding blank is really a proxy for what we are really interested in, which is the mean stimulus delivered to the neuron. In the following figure, we plot the gain (peak response/stimulus amplitude) as a function of the mean stimulus in the preceding 500ms.

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ph(1) = subplot(1,2,1); hold on 
axis square
ph(2) = subplot(1,2,2); hold on

s = 1; % when we start for the gain analysis
z = length(data(td).ORN); % where we end
example_history_length = 0.5;

clear x
x.response = data(td).ORN(s:z);
x.stimulus = data(td).PID(s:z);
x.time = data(td).time(s:z);
x.valve = data(td).Valve(s:z);

PulsePeakAnalysis(x,history_lengths,example_history_length,ph);


PrettyFig;

snapnow;
delete(gcf);	


%##       ##    ##    ##     ##  #######  ########  ######## ##          ########  ########     ##    
%##       ###   ##    ###   ### ##     ## ##     ## ##       ##          ##     ## ##     ##   ## #   
%##       ####  ##    #### #### ##     ## ##     ## ##       ##          ##     ## ##     ##  ##   #  
%##       ## ## ##    ## ### ## ##     ## ##     ## ######   ##          ########  ########  ##     # 
%##       ##  ####    ##     ## ##     ## ##     ## ##       ##          ##        ##        ######## 
%##       ##   ###    ##     ## ##     ## ##     ## ##       ##          ##        ##        ##     # 
%######## ##    ##    ##     ##  #######  ########  ######## ########    ##        ##        ##     ## 

%% Pulse Peak Analysis of LN Model
% In this section, we repeat the Peak-Pulse Analysis for the LN model prediction. 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ph(1) = subplot(1,2,1); hold on 
axis square
ph(2) = subplot(1,2,2); hold on

s = 1; % when we start for the gain analysis
z = length(data(td).ORN); % where we end
example_history_length = 0.5;

clear x
x.response = LNpred(s:z);
x.stimulus = data(td).PID(s:z);
x.time = data(td).time(s:z);
x.valve = data(td).Valve(s:z);

PulsePeakAnalysis(x,history_lengths,example_history_length,ph);

PrettyFig;


snapnow;
delete(gcf);



%% Version Info
% The file that generated this document is called:
disp(mfilename)

%%
% and its md5 hash is:
Opt.Input = 'file';
disp(DataHash(strcat(mfilename,'.m'),Opt))
