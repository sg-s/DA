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

% internal housekeeping: determine if being called by publish or not
calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end

history_lengths = (3*floor(1000*logspace(-2,1,30)/3))/1e3;
example_history_length = 0.135;


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
% The data looks like this. The following figure shows the valve state, the odor concentration, and the neuron response. The neuron is ab3A, and the odor presented is 1-octen-3-ol diluted to 3x $10^{-3}$ in Paraffin Oil. The correlation time in the valve position is 30ms. The stimulus measurement is de-trended with a quadratic term to remove what is supposed to be a drift in the PID. 

%%
% This data file is used for the following analysis:
td = 4;

disp(data(td).original_name)


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

if being_published
	snapnow;
	delete(fh);
end


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

if being_published
	snapnow;
	delete(gcf);
end

%%
% What does the prediction of this filter look like? The following figure shows the data (black) and the linear prediction (red) superimposed on it. 

figure('outerposition',[0 0 1100 400],'PaperUnits','points','PaperSize',[1100 400]); hold on
plot(data(td).time,data(td).ORN,'k');
plot(data(td).time,data(td).LinearFit,'r');
set(gca,'XLim',[min(data(td).time) max(data(td).time)])
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')

PrettyFig;

if being_published
	snapnow;
	delete(gcf);
end

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

if being_published
	snapnow;
	delete(gcf);
end


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
x = lsqcurvefit(@hill,[max(ydata) 2 2],xdata,ydata,[max(ydata)/2 2 1],[2*max(ydata) max(ydata) 10],fo);

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
plot(data(td).time,hill(x,data(td).LinearFit),'g')
set(gca,'XLim',[18 22])
xlabel('Time (s)')
ylabel('Firing rate (Hz)')

PrettyFig;

if being_published
	snapnow;
	delete(gcf);
end

% save this for later
LNpred = hill(x,data(td).LinearFit);

%%
% They look almost identical. The r-square of the LN prediction is: 

disp(rsquare(LNpred,data(td).ORN))

%%
% And the raw Euclidean Distance between the prediction and the data is: 

disp(Cost2(LNpred(205:end-33),data(td).ORN(205:end-33)))



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
	d.response = data(td).ORN;
	x_NLN = [.0068    2.447          130    4.76];
	[~,x_NLN] = FitNLNModel(d,x_NLN);
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
ms = 0; Ms = max(b);
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
set(gca,'XLim',[18 22])
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




%         ######      ###    #### ##    ##             ###    ##    ##    ###    ##       
%        ##    ##    ## ##    ##  ###   ##            ## ##   ###   ##   ## ##   ##       
%        ##         ##   ##   ##  ####  ##           ##   ##  ####  ##  ##   ##  ##       
%        ##   #### ##     ##  ##  ## ## ##          ##     ## ## ## ## ##     ## ##       
%        ##    ##  #########  ##  ##  ####          ######### ##  #### ######### ##       
%        ##    ##  ##     ##  ##  ##   ###          ##     ## ##   ### ##     ## ##       
%         ######   ##     ## #### ##    ##          ##     ## ##    ## ##     ## ######## 
      
%% Linear Model Gain Analysis
% In this section, we perform a gain analysis on the linear prediction first for an arbitrarily chosen history length of 120ms (panel on the left) and then for various history lengths (panel on the right). The history lengths where the slopes are significantly different (p<0.05) are indicated by dots. Significance is determined by bootstrapping the data 100 times. History lengths up to twice the autocorrelation length of the stimulus are investigated. 

%%
% We won't do this case because the fit is so bad. 

% data(td).LinearFit(data(td).LinearFit<0)=0;

% f1=figure('outerposition',[0 0 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
% ph(1) = subplot(2,1,1); hold on 
% ph(2) = subplot(2,1,2); hold on

% title(ph(1),strrep(data(td).original_name,'_','-'),'FontSize',20);

% f2=figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
% ph(3) = subplot(1,2,1); hold on 
% axis square
% ph(4) = subplot(1,2,2); hold on

% s = 300; % when we start for the gain analysis
% z = length(data(td).ORN) - 33; % where we end


% clear x
% x.response = data(td).ORN(s:z);
% x.prediction = data(td).LinearFit(s:z);
% x.stimulus = data(td).PID(s:z);
% x.time = data(td).time(s:z);
% x.filter_length = 201;


% if redo_bootstrap
% 	[p_K,l,h] = GainAnalysis4(x,history_lengths,example_history_length,ph);
% 	s=abs(l-h);
% 	s(p_K(1,:)>0.05)=NaN;
% 	[~,loc]=max(s);
% 	example_history_length_K = history_lengths(loc);
% else
% 	GainAnalysis4(x,history_lengths,example_history_length_K,ph,p_K);
% end

% xlabel(ph(3),'Linear Prediction (Hz)')
% set(ph(4),'XScale','log')

% if being_published
% 	snapnow;
% 	delete(f1);

% 	snapnow;
% 	delete(f2);
% end


%  ######      ###    #### ##    ##       ###    ##    ##    ###    ##          ##       ##    ## 
% ##    ##    ## ##    ##  ###   ##      ## ##   ###   ##   ## ##   ##          ##       ###   ## 
% ##         ##   ##   ##  ####  ##     ##   ##  ####  ##  ##   ##  ##          ##       ####  ## 
% ##   #### ##     ##  ##  ## ## ##    ##     ## ## ## ## ##     ## ##          ##       ## ## ## 
% ##    ##  #########  ##  ##  ####    ######### ##  #### ######### ##          ##       ##  #### 
% ##    ##  ##     ##  ##  ##   ###    ##     ## ##   ### ##     ## ##          ##       ##   ### 
%  ######   ##     ## #### ##    ##    ##     ## ##    ## ##     ## ########    ######## ##    ## 


%% LN Model Gain Analysis
% In the previous section, we saw that the simple linear model fails to account for this fast adaptation of the ORNs. Specifically, the gain of the neuron _w.r.t_ to the model is significantly different for times when the stimulus is high and the when the stimulus is low.

%%
% In this section, we want to know if the LN model (adding a non-linear function post-hoc) corrects for this mis-prediction of gain. Here, we repeat the gain analysis as in the previous section, but this time, using the LN prediction instead of the linear prediction. 

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

clear x
x.response = data(td).ORN(s:z);
x.prediction = LNpred(s:z);
x.stimulus = data(td).PID(s:z);
x.time = data(td).time(s:z);
x.filter_length = 201;

if redo_bootstrap
	[p_LN,l,h] = GainAnalysis4(x,history_lengths,example_history_length,ph);
	s=abs(l-h);
	s(p_LN(1,:)>0.05)=NaN;
	[~,loc]=max(s);
	example_history_length_LN = history_lengths(loc);

else
	GainAnalysis4(x,history_lengths,example_history_length_LN,ph,p_LN);
end

xlabel(ph(3),'LN Prediction (Hz)')
set(ph(4),'XScale','log')

if being_published
	snapnow;
	delete(f1);

	snapnow;
	delete(f2);
end

% ##    ## ##       ##    ##     ######      ###    #### ##    ##       ###    ##    ##    ###    ##       
% ###   ## ##       ###   ##    ##    ##    ## ##    ##  ###   ##      ## ##   ###   ##   ## ##   ##       
% ####  ## ##       ####  ##    ##         ##   ##   ##  ####  ##     ##   ##  ####  ##  ##   ##  ##       
% ## ## ## ##       ## ## ##    ##   #### ##     ##  ##  ## ## ##    ##     ## ## ## ## ##     ## ##       
% ##  #### ##       ##  ####    ##    ##  #########  ##  ##  ####    ######### ##  #### ######### ##       
% ##   ### ##       ##   ###    ##    ##  ##     ##  ##  ##   ###    ##     ## ##   ### ##     ## ##       
% ##    ## ######## ##    ##     ######   ##     ## #### ##    ##    ##     ## ##    ## ##     ## ######## 


%% NLN Model Gain Analysis
% In this section, we want to know if the NLN model shows the same features of fast gain control as we saw in the previous cases. Here, we repeat the gain analysis as in the previous section, but this time, using the NLN prediction instead of the linear prediction. 

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

clear x
x.response = data(td).ORN(s:z);
x.prediction = NLNFit(s:z);
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


%        ######      ###    #### ##    ##       ###    ##    ##    ###    ##       
%       ##    ##    ## ##    ##  ###   ##      ## ##   ###   ##   ## ##   ##       
%       ##         ##   ##   ##  ####  ##     ##   ##  ####  ##  ##   ##  ##       
%       ##   #### ##     ##  ##  ## ## ##    ##     ## ## ## ## ##     ## ##       
%       ##    ##  #########  ##  ##  ####    ######### ##  #### ######### ##       
%       ##    ##  ##     ##  ##  ##   ###    ##     ## ##   ### ##     ## ##       
%        ######   ##     ## #### ##    ##    ##     ## ##    ## ##     ## ######## 
      
%        ########  ######## ########    ###    #### ##        ######  
%       ##     ## ##          ##      ## ##    ##  ##       ##    ## 
%       ##     ## ##          ##     ##   ##   ##  ##       ##       
%       ##     ## ######      ##    ##     ##  ##  ##        ######  
%       ##     ## ##          ##    #########  ##  ##             ## 
%       ##     ## ##          ##    ##     ##  ##  ##       ##    ## 
%       ########  ########    ##    ##     ## #### ########  ######  


%% Understanding features of Linear Model Gain Analysis
% There are some features of the gain analysis that are interesting, especially when compared to the gain analysis of the NL model: 
% 
% * The slopes of the linear model gain analysis do not go to 1 even when the history lengths go to 0
% * The LN model does not behave this way: it behaves as expected, and the slopes go to zero as the history length go to 0
% * The green and the red curves are not mirror images of each other. Why not? What makes them different? 

%% 
% In the following figure, we use a DA model to generate synthetic responses to the actual stimulus and then fit a linear filter to this data, and then perform the gain analysis as above on this synthetic dataset. 

%% 
% The parameters of the DA model are chosen to be a best fit to the actual experimental data, so the DA model chosen best represents the actual neuron's response. The following figure shows the response of the neuron and the best-fit DA Model prediction. 



% if ~exist('p')
% 	d.stimulus = data(td).PID;
% 	d.response = data(td).ORN;
% 	p = FitDAModelToData(d,[25576 85 0.03 1.5 2 13 2 -151],[900 50 0 1e-1 2 1e-2 2 -200],[98000 5000 1 20 2 40 2 1],200);
% end



% figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
% f=DA_integrate2(data(td).PID,p);
% plot(data(td).time,data(td).ORN,'k')
% hold on
% plot(data(td).time,f,'r')
% xlabel('Time (s)')
% ylabel('Firing rate (Hz)')
% legend({'Data','DA Model'})
% set(gca,'XLim',[18 22])

% PrettyFig;

% if being_published
% 	snapnow;
% 	delete(gcf)
% end

% %%
% % The r-square of the DA Model fit is:
% disp(rsquare(data(td).ORN(300:end),f(300:end)))

% %% 
% % Now, we use the DA model as a "fake" neuron and use it to generate a synthetic ORN output, and then perform a linear gain analysis on this dataset. 

% % back out a filter
% s = 300; % when we start for the gain analysis
% z = length(data(td).ORN)-33; % where we end
% [K,~,filtertime] = FindBestFilter(data(td).PID,f);
% fp=mean(f(s:z))+convolve(data(td).time,data(td).PID,K,filtertime);



% f1=figure('outerposition',[0 0 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
% ph(1) = subplot(2,1,1); hold on 
% ph(2) = subplot(2,1,2); hold on

% title(ph(1),'Original Stimulus + DA Model Response','FontSize',20);

% f2=figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
% ph(3) = subplot(1,2,1); hold on 
% axis square
% ph(4) = subplot(1,2,2); hold on

% s = 300; % when we start for the gain analysis
% z = length(data(td).ORN) - 33; % where we end

% clear x
% x.response = f(s:z);
% x.prediction = fp(s:z);
% x.stimulus = data(td).PID(s:z);
% x.time = data(td).time(s:z);
% x.filter_length = 201;


% if redo_bootstrap
% 	ptemp3 = GainAnalysis4(x,history_lengths,example_history_length,ph);
% else
% 	GainAnalysis4(x,history_lengths,example_history_length,ph,ptemp3);
% end

% legend(ph(2),{'DA Response','Linear Prediction'})
% ylabel(ph(3),'DA Model Response (Hz)')
% xlabel(ph(3),'Linear Prediction (Hz)')
% set(ph(4),'XScale','log')

% if being_published
% 	snapnow;
% 	delete(f1);

% 	snapnow;
% 	delete(f2);
% end

% %%
% % So even the DA model can recapitulate the effect we see, where the slopes don't approach 1 as the history length goes to 0. So there's something about the way the linear prediction is made. 



%% Version Info
% The file that generated this document is called:
disp(mfilename)

%%
% and its md5 hash is:
Opt.Input = 'file';
disp(DataHash(strcat(mfilename,'.m'),Opt))