% Mahmut_Data_Analyis.m
% 
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

% specify which data to look at
load('/local-data/DA-paper/mahmut_data.mat')
td = 6;


roi1 = [20 30];
roi2 = [55 65];

% internal housekeeping: determine if being called by publish or not
calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end

example_history_length = 0.135;

%% 
% How do ORNs respond to non-Gaussian inputs? Real odor stimuli are characterised by long tails and non-gaussian statisitcs, with large whiffs of odor that occur in periods of relatively low signal. Such a stimulus has been generated here, and the responses of ORNs to these signals is analysed in this figure.





%                     ########     ###    ########    ###    
%                     ##     ##   ## ##      ##      ## ##   
%                     ##     ##  ##   ##     ##     ##   ##  
%                     ##     ## ##     ##    ##    ##     ## 
%                     ##     ## #########    ##    ######### 
%                     ##     ## ##     ##    ##    ##     ## 
%                     ########  ##     ##    ##    ##     ## 

%% Rough Overview of Data
% The following figure shows what the stimulus and the neuron's response looks like. 

% assuming baseline has already been removed...


figure('outerposition',[0 0 1000 800],'PaperUnits','points','PaperSize',[1000 800]); hold on
subplot(4,1,1), hold on
if isvector(data(td).PID)
	plot(data(td).time,data(td).PID,'k');
else
	plot(data(td).time,mean2(data(td).PID),'k');
end
ylabel('PID (V)')
titlestr = strkat(data(td).neuron,' response to ',data(td).odor);
title(titlestr)
set(gca,'XLim',[min(data(td).time) max(data(td).time)])
set(gca,'YLim',[-0.1 1.1*max(max(data(td).PID))])

subplot(4,1,2), hold on
if isvector(data(td).PID)
	plot(data(td).time,data(td).PID,'k');
else
	plot(data(td).time,mean2(data(td).PID),'k');
end
ylabel('PID (V)')
set(gca,'XLim',[min(data(td).time) max(data(td).time)])
set(gca,'YLim',[0 1.1*max(max(data(td).PID))])
a = min(nonzeros(data(td).PID(:)));
a = max([a 1e-6]);
z = 1.1*max((data(td).PID(:)));
set(gca,'YTick',logspace(log10(a),log10(z),4));
set(gca,'YScale','log','YMinorTick','on')

subplot(4,1,3), hold on
if isfield(data,'spiketimes') && isfield(data,'spiketimeB')
	if isempty(data(td).spiketimeB)
		raster2(data(td).spiketimes)
	else
		raster2(data(td).spiketimeB,data(td).spiketimes)
	end
end
set(gca,'XLim',[min(data(td).time) max(data(td).time)])

subplot(4,1,4), hold on
if isvector(data(td).ORN)
	plot(data(td).time,data(td).ORN,'k');
else
	plot(data(td).time,mean2(data(td).ORN),'k');
end
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
set(gca,'XLim',[min(data(td).time) max(data(td).time)])
set(gca,'YLim',[0 1.1*max(max(data(td).ORN))])

PrettyFig;

if being_published
	snapnow;
	delete(gcf);
end



%     ##     ##    ###    ########  ####    ###    ########  ##       #### ######## ##    ## 
%     ##     ##   ## ##   ##     ##  ##    ## ##   ##     ## ##        ##     ##     ##  ##  
%     ##     ##  ##   ##  ##     ##  ##   ##   ##  ##     ## ##        ##     ##      ####   
%     ##     ## ##     ## ########   ##  ##     ## ########  ##        ##     ##       ##    
%      ##   ##  ######### ##   ##    ##  ######### ##     ## ##        ##     ##       ##    
%       ## ##   ##     ## ##    ##   ##  ##     ## ##     ## ##        ##     ##       ##    
%        ###    ##     ## ##     ## #### ##     ## ########  ######## ####    ##       ##    

%% Trial to Trial variability
% How variable is the response and the stimulus from trial to trial? The following figure shows individual traces of the stimulus and the neuron's response for each trial. 


figure('outerposition',[0 0 1000 1000],'PaperUnits','points','PaperSize',[1000 1000]); hold on
set(gcf,'DefaultAxesColorOrder',jet(width(data(td).PID)))
subplot(3,2,1), hold on
plot(data(td).time,data(td).PID);
ylabel('PID (V)')
set(gca,'XLim',roi1)
set(gca,'YLim',[0 1.1*max(max(data(td).PID))])

subplot(3,2,2), hold on
plot(data(td).time,data(td).PID);
ylabel('PID (V)')
set(gca,'XLim',roi2)
set(gca,'YLim',[0 1.1*max(max(data(td).PID))])

if isfield(data,'spiketimes') && isfield(data,'spiketimeB')
	subplot(3,2,3), hold on
	raster2(data(td).spiketimes)
	set(gca,'XLim',roi1)


	subplot(3,2,4), hold on
	raster2(data(td).spiketimes)
	set(gca,'XLim',roi2)
end

subplot(3,2,5), hold on
plot(data(td).time,data(td).ORN);
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
set(gca,'XLim',roi1)
set(gca,'YLim',[0 1.1*max(max(data(td).ORN))])

subplot(3,2,6), hold on
plot(data(td).time,data(td).ORN);
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
set(gca,'XLim',roi2)
set(gca,'YLim',[0 1.1*max(max(data(td).ORN))])
PrettyFig;

if being_published
	snapnow;
	delete(gcf);
end



%%
% It looks like the amplitude of the signal is dropping every trial. It's still not clear why this is the case. From here onwards, the variability will be neglected, and all data from all the trials is averaged together. 
data(td).PID = mean2(data(td).PID);
data(td).ORN = mean2(data(td).ORN);

data(td).PID(data(td).PID<0) = 0;



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

if being_published
	snapnow;
	delete(gcf);
end




%     ##       #### ##    ## ########    ###    ########         ######## #### ######## 
%     ##        ##  ###   ## ##         ## ##   ##     ##        ##        ##     ##    
%     ##        ##  ####  ## ##        ##   ##  ##     ##        ##        ##     ##    
%     ##        ##  ## ## ## ######   ##     ## ########         ######    ##     ##    
%     ##        ##  ##  #### ##       ######### ##   ##          ##        ##     ##    
%     ##        ##  ##   ### ##       ##     ## ##    ##         ##        ##     ##    
%     ######## #### ##    ## ######## ##     ## ##     ##        ##       ####    ##    


%% Linear Fit
% The linear filter calculated and shown above is a "best" filter that simultaneously tries to maximise the predictive power of the model while suppressing high frequency components and having a unit gain. 

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

if being_published
	snapnow;
	delete(gcf);
end


%%
% The weird shape of the linear filter doesn't bode well for the quality of the linear prediction. The following figure shows the linear fit compared to the data. 

figure('outerposition',[0 0 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
subplot(2,2,1), hold on
if isvector(data(td).PID)
	plot(data(td).time,data(td).PID,'k');
else
	plot(data(td).time,mean2(data(td).PID),'k');
end
ylabel('PID (V)')
set(gca,'XLim',roi1)
set(gca,'YLim',[-0.1 1.1*max(max(data(td).PID))])

subplot(2,2,2), hold on
if isvector(data(td).PID)
	plot(data(td).time,data(td).PID,'k');
else
	plot(data(td).time,mean2(data(td).PID),'k');
end
ylabel('PID (V)')
set(gca,'XLim',roi2)
set(gca,'YLim',[-0.1 1.1*max(max(data(td).PID))])

subplot(2,2,3), hold on
if isvector(data(td).ORN)
	plot(data(td).time,data(td).ORN,'k');
else
	plot(data(td).time,mean2(data(td).ORN),'k');
end
plot(data(td).time,data(td).LinearFit,'r')
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
set(gca,'XLim',roi1)

subplot(2,2,4), hold on
if isvector(data(td).ORN)
	plot(data(td).time,data(td).ORN,'k');
else
	plot(data(td).time,mean2(data(td).ORN),'k');
end
plot(data(td).time,data(td).LinearFit,'r')
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
set(gca,'XLim',roi2)
legend Data LinearFit
PrettyFig;

if being_published
	snapnow;
	delete(fh);
end



%%
% The rsquare of the linear fit is:

if isvector(data(td).ORN)
	disp(rsquare(data(td).LinearFit,data(td).ORN))
else
	disp(rsquare(data(td).LinearFit,mean2(data(td).ORN)))
end


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
x = lsqcurvefit(@hill,[max(ydata) 2 2],xdata,ydata,[1 0 0],[2*max(ydata) max(ydata) 10],fo);
% save this for later
LNpred = hill(x,data(td).LinearFit);

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(xdata,hill(x,xdata),'k')
xlabel('Linear Prediction (Hz)')
ylabel('Nonlinearity Output (Hz)')

PrettyFig;

if being_published
	snapnow;
	delete(gcf);
end

%%
% That is almost linear. How does the output nonlinearity change the prediciton? In the following figure, the data is shown in black, the linear prediciton is shown in red, and the LN prediction is shown in green. 


figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(data(td).time,data(td).ORN,'k')
plot(data(td).time,data(td).LinearFit,'r')
plot(data(td).time,LNpred,'g')
set(gca,'XLim',[20 40])
xlabel('Time (s)')
ylabel('Firing rate (Hz)')
legend({'Data','Linear Fit','LN Fit'})

PrettyFig;

if being_published
	snapnow;
	delete(gcf);
end



%%
% They look almost identical. The r-square of the LN prediction is: 

disp(rsquare(LNpred,data(td).ORN))

%%
% And the raw Euclidean Distance between the prediction and the data is: 

disp(Cost2(LNpred(205:end-33),data(td).ORN(205:end-33)))


% taking log of stimulus
logPID = log(data(td).PID);
logPID= logPID - mean2(logPID(1:1000));
logPID(logPID<0) = 0;


%         ######      ###    #### ##    ##             ###    ##    ##    ###    ##       
%        ##    ##    ## ##    ##  ###   ##            ## ##   ###   ##   ## ##   ##       
%        ##         ##   ##   ##  ####  ##           ##   ##  ####  ##  ##   ##  ##       
%        ##   #### ##     ##  ##  ## ## ##          ##     ## ## ## ## ##     ## ##       
%        ##    ##  #########  ##  ##  ####          ######### ##  #### ######### ##       
%        ##    ##  ##     ##  ##  ##   ###          ##     ## ##   ### ##     ## ##       
%         ######   ##     ## #### ##    ##          ##     ## ##    ## ##     ## ######## 
      
%% log-NLN Model Gain Analysis
% In this section, we perform a gain analysis on the log NLN prediction first for an arbitrarily chosen history length of 120ms (panel on the left) and then for various history lengths (panel on the right). The history lengths where the slopes are significantly different (p<0.05) are indicated by dots. Significance is determined by bootstrapping the data 100 times.


% do it only when LNpred is non zero
CensoredLNPred = LNpred;
CensoredLNPred(LNpred<1) = NaN;
Censored_f = data(td).ORN;
Censored_f(LNpred<1) = NaN;

f1=figure('outerposition',[0 0 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
ph(1) = subplot(2,1,1); hold on 
ph(2) = subplot(2,1,2); hold on

f2=figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ph(3) = subplot(1,2,1); hold on 
axis square
ph(4) = subplot(1,2,2); hold on

s = 300; % when we start for the gain analysis
z = length(data(td).ORN) - 33; % where we end
history_lengths = (3*floor(1000*logspace(-2,1,30)/3))/1e3;

clear x
x.response = Censored_f(s:z);
x.prediction = CensoredLNPred(s:z);
x.stimulus = data(td).PID(s:z);
x.time = data(td).time(s:z);
x.filter_length = 201;


if redo_bootstrap
	[p_NLN,l,h] = GainAnalysis4(x,history_lengths,example_history_length,ph);
	s=abs(l-h);
	s(p_NLN(1,:)>0.05)=NaN;
	[~,plot_loc]=max(s);
	example_history_length_NLN = history_lengths(plot_loc);
else
	GainAnalysis4(x,history_lengths,example_history_length_NLN,ph,p_NLN);
end

xlabel(ph(3),'LN Prediction (Hz)')
set(ph(4),'XScale','log')

if being_published
	snapnow;
	delete(f1);

	snapnow;
	delete(f2);
end


return

% let's try to fit a DA model
clear x
x.response = data(td).ORN;
x.stimulus = data(td).PID;




%        ##    ## ##       ##    ##    ##     ##  #######  ########  ######## ##       
%        ###   ## ##       ###   ##    ###   ### ##     ## ##     ## ##       ##       
%        ####  ## ##       ####  ##    #### #### ##     ## ##     ## ##       ##       
%        ## ## ## ##       ## ## ##    ## ### ## ##     ## ##     ## ######   ##       
%        ##  #### ##       ##  ####    ##     ## ##     ## ##     ## ##       ##       
%        ##   ### ##       ##   ###    ##     ## ##     ## ##     ## ##       ##       
%        ##    ## ######## ##    ##    ##     ##  #######  ########  ######## ######## 


%% NLN Model
% In this section we fit the data with a NLN model, where there is a static nonlinear function both at the input and the output of the linear filter. The input and output non-linearities are two-parameter Hill functions. The linear filter is a non-parametric filter computed on-the-fly from the output of the first Hill function and the inverse of the second Hill function. All three parts of the model: the two functions and the filter -- are fitted simultaneously. 

%%
% The following figure shows the best fit NLN model.

if redo_bootstrap
	clear d
	d.stimulus = data(td).PID;
	d.stimulus(d.stimulus<0) = 0;
	d.response = data(td).ORN;
	x0 = [9.4453    0.9268   22.3047    1.5814];
	[NLNFit,x_NLN] = FitNLNModel(d,x0);
end

% solve the model for the best fit parameters
[NLNFit, K, a, b] = SolveNLNModel2(x_NLN,data(td).PID,data(td).ORN);

figure('outerposition',[0 0 1300 400],'PaperUnits','points','PaperSize',[1300 400]); hold on
subplot(1,3,1), hold on
ms = min(data(td).PID); Ms = max(data(td).PID);
xx = ms:(Ms-ms)/100:Ms;
y = hill2(x_NLN(1:2),xx);
plot(xx,y)
xlabel('Stimulus (V)')
ylabel('Input to filter (V)')

subplot(1,3,2), hold on
t = 1:300;
tt = 3e-3:3e-3:3e-3*length(K);
plot(tt,K);
set(gca,'XLim',[min(tt) max(tt)])
xlabel('Filter Lag (s)')
ylabel('Filter Amplitude (Hz)')

b(isinf(b)) = [];
subplot(1,3,3), hold on
ms = 0; Ms = 300;
xx = ms:(Ms-ms)/100:Ms;
y = hill2(x_NLN(3:4),xx);
plot(xx,y)
xlabel('Filter output (Hz)')
ylabel('Model Output (Hz)')

PrettyFig;

if being_published
	snapnow;
	delete(gcf);
end

%%
% The following figure compares the NLN Model fit to the data:

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(data(td).time,data(td).ORN,'k')
plot(data(td).time,NLNFit,'r')
set(gca,'XLim',[20 40])
xlabel('Time (s)')
ylabel('Firing rate (Hz)')
legend({'Data','NLN Fit'})

PrettyFig;

if being_published
	snapnow;
	delete(gcf);
end

%%
%  The r-square of the NLN prediction is: 

disp(rsquare(NLNFit,data(td).ORN))

%%
% And the raw Euclidean Distance between the prediction and the data is: 

disp(Cost2(NLNFit(205:end-33),data(td).ORN(205:end-33)))



%     ##        #######   ######       ######  ######## ##    ##  ######  #### ##    ##  ######   
%     ##       ##     ## ##    ##     ##    ## ##       ###   ## ##    ##  ##  ###   ## ##    ##  
%     ##       ##     ## ##           ##       ##       ####  ## ##        ##  ####  ## ##        
%     ##       ##     ## ##   ####     ######  ######   ## ## ##  ######   ##  ## ## ## ##   #### 
%     ##       ##     ## ##    ##           ## ##       ##  ####       ##  ##  ##  #### ##    ##  
%     ##       ##     ## ##    ##     ##    ## ##       ##   ### ##    ##  ##  ##   ### ##    ##  
%     ########  #######   ######       ######  ######## ##    ##  ######  #### ##    ##  ######   

%% Logarithmic Sensing + NLN Model
% Even the NLN model seems to do a very poor job of predicting the response of the neuron. In this dataset, the stimulus varies wildly (by design), and we know that neurons respond logarithmically to the stimulus (the dose-response curve is log-linear). So perhaps we can predict the response better if we take the _logarithm_ of the stimulus, instead of the stimulus itself. 




%%
% The following figure shows the best fit NLN model.

if redo_bootstrap
	clear d
	d.stimulus = logPID;
	d.stimulus(d.stimulus<0) = 0;
	d.response = data(td).ORN;
	x0 = [8.3175    6.0375   48.4250    2.1891];
	[~,x_logNLN] = FitNLNModel(d,x0);
end

% solve the model for the best fit parameters
[logNLNFit, K, a, b] = SolveNLNModel2(x_logNLN,logPID,data(td).ORN);

figure('outerposition',[0 0 1300 400],'PaperUnits','points','PaperSize',[1300 400]); hold on
subplot(1,3,1), hold on
ms = min(logPID); Ms = max(logPID);
xx = ms:(Ms-ms)/100:Ms;
y = hill2(x_logNLN(1:2),xx);
plot(xx,y)
xlabel('Stimulus (V)')
ylabel('Input to filter (V)')

subplot(1,3,2), hold on
t = 1:300;
tt = 3e-3:3e-3:3e-3*length(K);
plot(tt,K);
set(gca,'XLim',[min(tt) max(tt)])
xlabel('Filter Lag (s)')
ylabel('Filter Amplitude (Hz)')

b(isinf(b)) = [];
subplot(1,3,3), hold on
ms = 0; Ms = 300;
xx = ms:(Ms-ms)/100:Ms;
y = hill2(x_logNLN(3:4),xx);
plot(xx,y)
xlabel('Filter output (Hz)')
ylabel('Model Output (Hz)')

PrettyFig;

if being_published
	snapnow;
	delete(gcf);
end

%%
% The following figure compares the log-NLN Model fit to the data:

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(data(td).time,data(td).ORN,'k')
plot(data(td).time,logNLNFit,'r')
set(gca,'XLim',[20 40])
xlabel('Time (s)')
ylabel('Firing rate (Hz)')
legend({'Data','log-NLN Fit'})

PrettyFig;

if being_published
	snapnow;
	delete(gcf);
end

%%
%  The r-square of the NLN prediction is: 

disp(rsquare(logNLNFit,data(td).ORN))

%%
% And the raw Euclidean Distance between the prediction and the data is: 

disp(Cost2(logNLNFit(205:end-33),data(td).ORN(205:end-33)))


%    ######## #### ##     ## ######## ########     ######## #### ##       ######## ######## ########  
%    ##        ##   ##   ##  ##       ##     ##    ##        ##  ##          ##    ##       ##     ## 
%    ##        ##    ## ##   ##       ##     ##    ##        ##  ##          ##    ##       ##     ## 
%    ######    ##     ###    ######   ##     ##    ######    ##  ##          ##    ######   ########  
%    ##        ##    ## ##   ##       ##     ##    ##        ##  ##          ##    ##       ##   ##   
%    ##        ##   ##   ##  ##       ##     ##    ##        ##  ##          ##    ##       ##    ##  
%    ##       #### ##     ## ######## ########     ##       #### ########    ##    ######## ##     ## 

%% Fixing the filters
% This data set is not ideal for extracting a linear filter from. Specifically, we need some sort of "white noise" like stimuli to correctly estimate a filters. So we are now going to use an "average" filter from the flickering stimulus experiment and fix that in the NLN model.  



load('HowGeneralIsGainAdaptation.mat','Filters','filtertime');
K = mean2(Filters(:,2:end));
K(filtertime<0) = [];
K = K/max(K);
filtertime(filtertime<0) = [];

if redo_bootstrap
	clear d
	d.stimulus = logPID;
	d.stimulus(d.stimulus<0) = 0;
	d.response = data(td).ORN;
	d.K = K;
	x0 = [9.5146    2.7273   43.9182    1.6673];
	[~,x_fixedK] = FitNLNModel(d,x0);
end


% solve the model for the best fit parameters
[fixedKFit, ~, a, b] = SolveNLNModel2(x_fixedK,logPID,data(td).ORN,K);

figure('outerposition',[0 0 1300 400],'PaperUnits','points','PaperSize',[1300 400]); hold on
subplot(1,3,1), hold on
ms = min(logPID); Ms = max(logPID);
xx = ms:(Ms-ms)/100:Ms;
y = hill2(x_fixedK(1:2),xx);
plot(xx,y)
xlabel('Stimulus (V)')
ylabel('Input to filter (V)')

subplot(1,3,2), hold on
t = 1:300;
tt = 3e-3:3e-3:3e-3*length(K);
plot(tt,K);
set(gca,'XLim',[min(tt) max(tt)])
xlabel('Filter Lag (s)')
ylabel('Filter Amplitude (Hz)')

b(isinf(b)) = [];
subplot(1,3,3), hold on
ms = 0; Ms = 300;
xx = ms:(Ms-ms)/100:Ms;
y = hill2(x_fixedK(3:4),xx);
plot(xx,y)
xlabel('Filter output (Hz)')
ylabel('Model Output (Hz)')

PrettyFig;

if being_published
	snapnow;
	delete(gcf);
end

%%
% The following figure compares the log-NLN Model fit to the data:

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(data(td).time,data(td).ORN,'k')
plot(data(td).time,fixedKFit,'r')
set(gca,'XLim',[20 40])
xlabel('Time (s)')
ylabel('Firing rate (Hz)')
legend({'Data','fixed-K Fit'})

PrettyFig;

if being_published
	snapnow;
	delete(gcf);
end

%%
%  The r-square of the NLN prediction is: 

disp(rsquare(fixedKFit,data(td).ORN))

%%
% And the raw Euclidean Distance between the prediction and the data is: 

disp(Cost2(fixedKFit(205:end-33),data(td).ORN(205:end-33)))

return


%%
% It is clear from the diagram that our method of gain analysis is all wrong: we need to constrain our analysis to the whiffs of odor, and perhaps compare large whiffs to small whiffs. 


%    ##      ## ##     ## #### ######## ########          ###    ##    ##    ###    ##       
%    ##  ##  ## ##     ##  ##  ##       ##               ## ##   ###   ##   ## ##   ##       
%    ##  ##  ## ##     ##  ##  ##       ##              ##   ##  ####  ##  ##   ##  ##       
%    ##  ##  ## #########  ##  ######   ######         ##     ## ## ## ## ##     ## ##       
%    ##  ##  ## ##     ##  ##  ##       ##             ######### ##  #### ######### ##       
%    ##  ##  ## ##     ##  ##  ##       ##             ##     ## ##   ### ##     ## ##       
%     ###  ###  ##     ## #### ##       ##             ##     ## ##    ## ##     ## ######## 

%% Whiff-based Analysis
% In this section, we consider the gain of the neuron on a whiff-by-whiff basis. To do so, we must identify the whiffs accurately. In the following figure, the whiff peaks (black), the whiff starts (green) and the whiff ends (red) are indicated. These were calculated automatically. 


p = data(td).PID;

figure('outerposition',[0 0 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
plot(data(td).time, p,'k'), hold on

% find peaks

m = mean(p(1:2000));
s = std(p(1:2000));
[pp,loc] = findpeaks(p,'MinPeakHeight',3*s+m,'MinPeakDistance',30);

% exclude false peaks created by stupid algo
rm_these = [];
sf = 1;
w = 20;
for i =1:length(loc)
	if mean(p(loc(i)-w:loc(i)-1)) < sf*pp(i)  && mean(p(loc(i)+1:loc(i)+w)) < sf*pp(i)
	else
		rm_these = [i rm_these];
	end

end
pp(rm_these)= [];
loc(rm_these) = [];
scatter(data(td).time(loc),p(loc),'b','filled');


% find whiff starts and stops
whiff_ons = 0*pp;
whiff_offs = 0*pp;
temp  =[]; temp2 = [];
for i = 1:length(loc)
	if i > 1
		a = loc(i-1);
	else
		a=1;
	end
	temp =  find(p(a:loc(i))<m+3*s,1,'last');

	if isempty(temp)
		% find the minimum
		[~,temp2]=min(p(a:loc(i)));
		whiff_ons(i) = temp2+a;
	else
		% all OKi
		whiff_ons(i)  = temp+a;
	end


end
for i = 1:length(loc)
	if i ==length(loc)
		z = length(p);
	else
		z = whiff_ons(i+1);
	end
	temp =  find(p(loc(i):z)<m+3*s,1,'first');
	if isempty(temp)
		% find the minimum
		[~,temp2]=min(p(loc(i):z));
		whiff_offs(i) = temp2+loc(i);
	else
		% all OKi
		whiff_offs(i)  = temp+loc(i);
	end

end
clear temp temp2
scatter(data(td).time(whiff_ons),p(whiff_ons),'g','filled');
scatter(data(td).time(whiff_offs),p(whiff_offs),'r','filled');

xlabel('Time (s)')
ylabel('Stimulus (V)')
title('Whiff identification')
set(gca,'XLim',[57 66],'YLim',[-.1 max(p)+1])

PrettyFig;

if being_published
	snapnow;
	delete(gcf);
end

% assemble the whiff data
whiff = struct;
for i = 1:length(whiff_ons)
	whiff(i).t_on = whiff_ons(i);
	whiff(i).t_off = whiff_offs(i);
	whiff(i).t_peak = loc(i);
	whiff(i).f = data(td).ORN(whiff_ons(i):whiff_offs(i));
	whiff(i).fp = logNLNFit(whiff_ons(i):whiff_offs(i));
	% fit a line 
	[fitmetrics, gof] = fit(whiff(i).fp,whiff(i).f,'Poly1');
	whiff(i).slope = fitmetrics.p1;
	whiff(i).gof = gof.rsquare;
	whiff(i).r = max(whiff(i).f)/max(whiff(i).fp);
end


% now only do gain analysis for the whiffs
WhiffORN = NaN*data(td).ORN;
WhiffORN_Pred = NaN*data(td).ORN;

for i = 1:length(whiff)
	WhiffORN(whiff(i).t_on:whiff(i).t_off) = whiff(i).f; 
	WhiffORN_Pred(whiff(i).t_on:whiff(i).t_off) = whiff(i).fp; 
end

f1=figure('outerposition',[0 0 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
ph(1) = subplot(2,1,1); hold on 
ph(2) = subplot(2,1,2); hold on

f2=figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ph(3) = subplot(1,2,1); hold on 
axis square
ph(4) = subplot(1,2,2); hold on

s = 300; % when we start for the gain analysis
z = length(data(td).ORN) - 33; % where we end
history_lengths = (3*floor(1000*logspace(-2,1,30)/3))/1e3;

clear x
x.response = WhiffORN(s:z);
x.prediction = WhiffORN_Pred(s:z);
x.stimulus = data(td).PID(s:z);
x.time = data(td).time(s:z);
x.filter_length = 201;


if redo_bootstrap
	[p_NLN,l,h] = GainAnalysis4(x,history_lengths,example_history_length,ph);
	s=abs(l-h);
	s(p_NLN(1,:)>0.05)=NaN;
	[~,loc]=max(s);
	example_history_length_NLN = history_lengths(loc);
else
	GainAnalysis4(x,history_lengths,example_history_length_NLN,ph,p_NLN);
end

xlabel(ph(3),'NLN Prediction (Hz)')
set(ph(4),'XScale','log')

if being_published
	snapnow;
	delete(f1);

	snapnow;
	delete(f2);
end













whiff([whiff.gof]<.8) = [];

g = NaN*history_lengths;
hl = floor(history_lengths/3e-3);
shat = ComputeSmoothedStimulus(data(td).PID,hl);
for i = 1:length(g)
	x = shat(i,[whiff.t_peak]); x = x(:);
	y = [whiff.slope]; y = y(:);
	y(isnan(x)) = []; x(isnan(x)) = [];
	[~,temp] = fit(x,y,'Poly1');
	g(i) = temp.rsquare;
end

return





%%
% Now we break up the trace so that we can perform a whiff-by-whiff analysis of the stimulus and the response, for each trial. The following figure shows how the whiff amplitude and response amplitude drop as a function of trial number.


if whiff_anal

	[whiff_ons,whiff_offs] = ComputeOnsOffs(hs);

	whiff_durations = whiff_offs - whiff_ons;
	whiff_ons(whiff_durations <20) = [];
	whiff_offs(whiff_durations <20) = [];
	whiff_durations = whiff_offs - whiff_ons;

	whiff_stim_max  = NaN(width(data(td).PID),length(whiff_ons));
	whiff_resp_max  = NaN(width(data(td).PID),length(whiff_ons));


	for j = 1:width(data(td).PID)
		for i = 1:length(whiff_ons)
			whiff_stim_max(j,i) = max(data(td).PID(whiff_ons(i):whiff_offs(i),j));
			whiff_resp_max(j,i) = max(data(td).ORN(whiff_ons(i):whiff_offs(i),j));
		end
		clear i

	end
	clear j

	for i = 1:length(whiff_ons)
		% normalise
		whiff_stim_max(:,i) = whiff_stim_max(:,i)/(whiff_stim_max(1,i));
		whiff_resp_max(:,i) = whiff_resp_max(:,i)/(whiff_resp_max(1,i));
	end

	figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
	set(gcf,'DefaultAxesColorOrder',jet(length(whiff_ons)))
	subplot(1,2,1), hold on
	plot(whiff_stim_max)
	xlabel('Trial Number')
	ylabel('Whiff Peak (norm)')
	set(gca,'YLim',[0 2])

	subplot(1,2,2), hold on
	plot(whiff_resp_max)
	xlabel('Trial Number')
	ylabel('Max Whiff Response (norm)')
	set(gca,'YLim',[0 2])

	PrettyFig;
	snapnow;
	delete(gcf);

end







return




%      ########     ###             ##     ##  #######  ########  ######## ##       
%      ##     ##   ## ##            ###   ### ##     ## ##     ## ##       ##       
%      ##     ##  ##   ##           #### #### ##     ## ##     ## ##       ##       
%      ##     ## ##     ##          ## ### ## ##     ## ##     ## ######   ##       
%      ##     ## #########          ##     ## ##     ## ##     ## ##       ##       
%      ##     ## ##     ##          ##     ## ##     ## ##     ## ##       ##       
%      ########  ##     ##          ##     ##  #######  ########  ######## ######## 



%% Fitting a DA Model
% Can a DA Model explain the responses of this neuron in this dataset? The following figure shows the ORN firing rates and the best-fit DA Model. 

if isvector(data(td).PID)
	data(td).DAFit=DA_integrate2(data(td).PID,data(td).DAFitParam);
else
	data(td).DAFit=DA_integrate2(mean2(data(td).PID),data(td).DAFitParam);
end


figure('outerposition',[0 0 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
subplot(2,2,1), hold on
if isvector(data(td).PID)
	plot(data(td).time,data(td).PID,'k');
else
	plot(data(td).time,mean2(data(td).PID),'k');
end
ylabel('PID (V)')
set(gca,'XLim',[25 35])
set(gca,'YLim',[0 1.1*max(max(data(td).PID))])

subplot(2,2,2), hold on
if isvector(data(td).PID)
	plot(data(td).time,data(td).PID,'k');
else
	plot(data(td).time,mean2(data(td).PID),'k');
end
ylabel('PID (V)')
set(gca,'XLim',[38 48])
set(gca,'YLim',[0 1.1*max(max(data(td).PID))])

subplot(2,2,3), hold on
if isvector(data(td).ORN)
	plot(data(td).time,data(td).ORN,'k');
else
	plot(data(td).time,mean2(data(td).ORN),'k');
end
plot(data(td).time,data(td).DAFit,'r')
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
set(gca,'XLim',[25 35])

subplot(2,2,4), hold on
if isvector(data(td).PID)
	plot(data(td).time,data(td).ORN,'k');
else
	plot(data(td).time,mean2(data(td).ORN),'k');
end
plot(data(td).time,data(td).DAFit,'r')
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
set(gca,'XLim',[38 48])
legend Data DAFit
PrettyFig;



snapnow;
delete(gcf);


%%
% The r-square of the DA model fit is:

if isvector(data(td).ORN)
	disp(rsquare(data(td).DAFit,data(td).ORN))
else
	disp(rsquare(data(td).DAFit,mean2(data(td).ORN)))
end

%%
% and the parameters of the model are:

disp(data(td).DAFitParam)


%% Gain Analysis: Comparison to DA Model
% Does the DA model also show systematic variation of gain for high and low stimuli? 

figure('outerposition',[0 0 1000 400],'PaperUnits','points','PaperSize',[1300 400]); hold on
clear ph
ph(3)=subplot(1,2,1); hold on; 	axis square
ph(4)=subplot(1,2,2); hold on;	axis square
s = 1; % when we start for the gain analysis
z = length(data(td).ORN); % where we end
step_size = 2*act/10;
if step_size < 0.03
	step_size= 0.03;
else
	step_size = 0.03*floor(step_size/0.03);
end
history_lengths = 0:step_size:2*act;
example_history_length = history_lengths(3);

clear x
if isvector(data(td).ORN)
	x.response = data(td).ORN(s:z);
else
	x.response = mean(data(td).ORN(s:z,:),2);
end
x.prediction = data(td).DAFit(s:z);
if isvector(data(td).PID)
	x.stimulus = data(td).PID(s:z);
else
	x.stimulus = mean(data(td).PID(s:z,:),2);
end
x.time = data(td).time(s:z);
x.filter_length = 201;

redo_bootstrap = 0;
if redo_bootstrap
	data(td).LinearFit_p = GainAnalysis3(x,history_lengths,example_history_length,ph);
else
	GainAnalysis3(x,history_lengths,example_history_length,ph,NaN*history_lengths);
end
clear x


snapnow;
delete(gcf);

%%
% Again, we repeat the gain analysis but restrict the analysis only to the whiffs, ignoring the blanks inbetween. 


if whiff_anal

	figure('outerposition',[0 0 1000 400],'PaperUnits','points','PaperSize',[1300 400]); hold on
	clear ph
	ph(3)=subplot(1,2,1); hold on; 	axis square
	ph(4)=subplot(1,2,2); hold on;	axis square
	s = 1; % when we start for the gain analysis
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
	if isvector(data(td).ORN)
		x.response = data(td).ORN(s:z);
	else
		x.response = mean(data(td).ORN(s:z,:),2);
	end
	x.prediction = data(td).DAFit(s:z);
	if isvector(data(td).PID)
		x.stimulus = data(td).PID(s:z);
	else
		x.stimulus = mean(data(td).PID(s:z,:),2);
	end

	% restrict to whiffs
	x.stimulus(~hs) = NaN;

	x.time = data(td).time(s:z);
	x.filter_length = 201;

	redo_bootstrap = 0;
	if redo_bootstrap
		data(td).LinearFit_p = GainAnalysis3(x,history_lengths,example_history_length,ph);
	else
		GainAnalysis3(x,history_lengths,example_history_length,ph,NaN*history_lengths);
	end
	clear x



	snapnow;
	delete(gcf);



end




%% Version Info
% The file that generated this document is called:
disp(mfilename)

%%
% and its md5 hash is:
Opt.Input = 'file';
disp(DataHash(strcat(mfilename,'.m'),Opt))

	