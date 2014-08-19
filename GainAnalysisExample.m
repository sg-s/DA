% GainAnalysisExample.m
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


% some parameters
font_size = 20;
marker_size = 10;
marker_size2 = 20;

redo_bootstrap = 0;

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
% The data looks like this. The following figure shows the valve state, the odor concentration, and the neuron response. The neuron is ab3A, and the odor presented is 1-octen-3-ol diluted to 3x $10^{-3}$ in Paraffin Oil. The correlation time in the valve position is 30ms. The stimulus measurement is de-trended with a quadratic term to remove what is supposed to be a drift in the PID. 

%%
% This data file is used for the following analysis:
td = 7;

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
% Once again, we need to correct the linear prediction by adding back the mean of the response: 

data(td).LinearFit = data(td).LinearFit + mean(data(td).ORN);

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
% That is almost linear. How does the output nonlinearity change the prediciton? In the following figure, the data is shown in black, the linear prediciton is shown in red, and the LN prediction is shown in green. 


figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(data(td).time,data(td).ORN,'k')
plot(data(td).time,data(td).LinearFit,'r')
plot(data(td).time,hill(x,data(td).LinearFit),'g')
set(gca,'XLim',[18 22])
xlabel('Time (s)')
ylabel('Firing rate (Hz)')

PrettyFig;

snapnow;
delete(gcf);

% save this for later
LNpred = hill(x,data(td).LinearFit);

%%
% They look almost identical. The r-square of the LN prediction is: 

disp(rsquare(hill(x,data(td).LinearFit),data(td).ORN))

%%
% And the raw Euclidian Distance between the prediction and the data is: 

disp(Cost2(hill(x,data(td).LinearFit(205:end-33)),data(td).ORN(205:end-33)))


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



%% LN Model Gain Analysis
% In the previous section, we saw that the simple linear model fails to account for this fast adaptation of the ORNs. Specifically, the gain of the neuron w.r.t to the model is significantly differnet for times when the stimulus is high and the when the stimulus is low.

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


