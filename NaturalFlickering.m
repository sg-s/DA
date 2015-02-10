% Mahmut_Data_Analyis.m
% this is a complete rewrite on  4:06 , 10 February 2015. 
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end

load('/local-data/DA-paper/natural-flickering/mahmut-raw/2014_07_11_EA_natflick_non_period_CFM_1_ab3_1_1_all.mat')
PID = data(2).PID;
time = 1e-4*(1:length(PID));
all_spikes = spikes(2).A;
B_spikes = spikes(2).B;


% A spikes --> firing rate
hash = DataHash(full(all_spikes));
cached_data = cache(hash);
if isempty(cached_data)
	fA = spiketimes2f(all_spikes,time);
	cache(hash,fA);
else
	fA = cached_data;
end

% B spikes --> firing rate
hash = DataHash(full(B_spikes));
cached_data = cache(hash);
if isempty(cached_data)
	fB = spiketimes2f(B_spikes,time);
	cache(hash,fB);
else
	fB = cached_data;
end

tA = 1e-3*(1:length(fA));
PID2 = fA;
for i = 1:width(PID2)
	PID2(:,i) = interp1(time,PID(i,:),tA);
end
PID = PID2; clear PID2
% some minor cleaning up
PID(end,:) = PID(end-1,:); 



figure('outerposition',[0 0 1400 1000],'PaperUnits','points','PaperSize',[1400 1000]); hold on
subplot(2,9,10:14), hold on
plot(tA,mean2(fA),'k')
set(gca,'XLim',[0 70])
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')

subplot(2,9,15:16), hold on
[r2,s] = rsquare(fA);
imagescnan(r2)
caxis([0 1])
axis image
colorbar
axis off
title(strcat('mean r^2 = ',oval(mean(r2(~isnan(r2))),2)))

subplot(2,9,17:18), hold on
imagescnan(s)
colorbar
axis image
axis off
title(strcat('mean slope = ',oval(mean(s(~isnan(s))),2)))


subplot(2,9,1:5), hold on
plot(tA,mean2(PID),'k')
set(gca,'YScale','log')
set(gca,'XLim',[0 70])
xlabel('Time (s)')
ylabel('Odor Concentration (V)')

subplot(2,9,6:7), hold on
[r2,s] = rsquare(PID);
imagescnan(r2)
caxis([0 1])
colorbar
axis image
axis off
title(strcat('mean r^2 = ',oval(mean(r2(~isnan(r2))),2)))

subplot(2,9,8:9), hold on
imagescnan(s)
colorbar
axis image
axis off
title(strcat('mean slope = ',oval(mean(s(~isnan(s))),2)))

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%%
% For completeness, here is a comparison of the A neuron's response to that of the B neuron's. 


figure('outerposition',[0 0 1400 1000],'PaperUnits','points','PaperSize',[1400 1000]); hold on
subplot(2,8,1:6), hold on
plot(tA,mean2(fA),'k')
set(gca,'XLim',[0 70])
xlabel('Time (s)')
ylabel('Firing Rate (A) (Hz)')
set(gca,'YLim',[0 120])

subplot(2,8,7:8), hold on
hash = DataHash(fA);
cached_data = cache(hash);
if isempty(cached_data)
	r2 = rsquare(fA);
	cache(hash,r2);
else
	r2 = cached_data;
end
imagescnan(r2)
caxis([0 1])
colorbar
axis image
axis off
title(strcat('mean r^2 = ',oval(mean(r2(~isnan(r2))),2)))


subplot(2,8,9:14), hold on
plot(tA,mean2(fB),'k')
set(gca,'XLim',[0 70])
xlabel('Time (s)')
ylabel('Firing Rate (B) (Hz)')
set(gca,'YLim',[0 120])

subplot(2,8,15:16), hold on
hash = DataHash(fB);
cached_data = cache(hash);
if isempty(cached_data)
	r2 = rsquare(fB);
	cache(hash,r2);
else
	r2 = cached_data;
end
imagescnan(r2)
caxis([0 1])
colorbar
axis image
axis off
title(strcat('mean r^2 = ',oval(mean(r2(~isnan(r2))),2)))


PrettyFig();

if being_published
	snapnow
	delete(gcf)
end



%%
% In this document, we generate odor stimuli flickers over a large range, like the "natural" stimuli, but never goes to zero, so that the neuron should never silence (allowing us to accurately follow its response, and reasonably estimate instantaneous gain). 

%%
% The first figure shows the odor stimulus (ethyl acetate) presented to two ab3A neurons. Each neuron was recorded from ten times. The colormaps on the right show the coefficient of determination between pairwise trials. 


%%
% The following figure shows the stimulus distribution and the autocorrelation functions of the stimulus and the response. 


figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1), hold on
for i = 1:width(PID)
	[y,x] = histcounts(PID(:,i),300);x(1) = [];
	plot(x,y)
end

set(gca,'YScale','log','XScale','log')
xlabel('PID (V)')
ylabel('Count')
title('Stimulus Histogram')

subplot(1,3,2), hold on
a = [];
for i = 1:width(PID)
	[~,~,c,temp]=FindCorrelationTime(PID(:,i));
	a = [a temp];
	l=plot(1:3e3,c(1:3e3));
end

xlabel('Lag (ms)')
ylabel('Autocorrelation')
legend(l,strcat('\tau=',oval(a),'ms'))
set(gca,'XScale','log')
title('Stimulus')

subplot(1,3,3), hold on
a = [];
for i = 1:width(fA)
	[~,~,c,temp]=FindCorrelationTime(fA(:,i));
	a = [a temp];
	l=plot(1:3e3,c(1:3e3));
end

title('Response')
xlabel('Lag (ms)')
ylabel('Autocorrelation')
legend(l,strcat('\tau=',oval(a),'ms'))
set(gca,'XScale','log')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end



%     ##       #### ##    ## ########    ###    ########         ######## #### ######## 
%     ##        ##  ###   ## ##         ## ##   ##     ##        ##        ##     ##    
%     ##        ##  ####  ## ##        ##   ##  ##     ##        ##        ##     ##    
%     ##        ##  ## ## ## ######   ##     ## ########         ######    ##     ##    
%     ##        ##  ##  #### ##       ######### ##   ##          ##        ##     ##    
%     ##        ##  ##   ### ##       ##     ## ##    ##         ##        ##     ##    
%     ######## #### ##    ## ######## ##     ## ##     ##        ##       ####    ##    


%% Linear Fit
% In this section we fit a linear kernel, first to the stimulus, and then to the log of the stimulus using standard methods. 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[900 500]); hold on
subplot(1,2,1), hold on
plot([-.1 1],[0 0 ],'k--')
[K, ~, filtertime_full] = FindBestFilter(mean2(PID),mean2(fA),[],'regmax=1;','regmin=1;','filter_length=1999;','offset=500;');
filtertime_full = filtertime_full*mean(diff(tA));
filtertime = 1e-3*(-200:900);
K = interp1(filtertime_full,K,filtertime);

plot(filtertime,K,'r')
ylabel('Filter Amplitude')
xlabel('Filter Lag (s)')
title('Filter from stimulus')

subplot(1,2,2), hold on
plot([-.1 1],[0 0 ],'k--')
[Klog, ~, filtertime_full] = FindBestFilter(log(mean2(PID)),mean2(fA),[],'regmax=1;','regmin=1;','filter_length=1999;','offset=500;');
filtertime_full = filtertime_full*mean(diff(tA));
filtertime = 1e-3*(-200:900);
Klog = interp1(filtertime_full,Klog,filtertime);

plot(filtertime,Klog,'r')
ylabel('Filter Amplitude')
xlabel('Filter Lag (s)')
title('Filter from log stimulus')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end


%%
% We will use this filter to make a prediction of the response:

% convolve with filter to make prediction
fp = convolve(tA,mean2(PID),K,filtertime);
fp_log = convolve(tA,log(mean2(PID)),Klog,filtertime);

% correct for trivial scaling
R = mean2(fA);
temp =fit(fp(~(isnan(fp) | isnan(R))),R(~(isnan(fp) | isnan(R))),'poly1');
fp = fp*temp.p1;
fp = fp+temp.p2;
temp =fit(fp_log(~(isnan(fp_log) | isnan(R))),R(~(isnan(fp_log) | isnan(R))),'poly1');
fp_log = fp_log*temp.p1;
fp_log = fp_log+temp.p2;


figure('outerposition',[0 0 1500 800],'PaperUnits','points','PaperSize',[1500 800]); hold on
subplot(2,4,1:3), hold on
plot(tA,mean2(fA),'k')
l=plot(tA,fp,'r');
r2 = rsquare(fp,mean2(fA));
legend(l,strcat('r^2=',oval(r2,2)))
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')

subplot(2,4,4), hold on
plot(filtertime,K,'r')
xlabel('Filter Lag (s)')
ylabel('Filter (norm.)')
title('Filter from stimulus')

subplot(2,4,5:7), hold on
plot(tA,mean2(fA),'k')
l=plot(tA,fp_log,'r');
r2 = rsquare(fp_log,mean2(fA));
legend(l,strcat('r^2=',oval(r2,2)))
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')

subplot(2,4,8), hold on
plot(filtertime,K,'r')
xlabel('Filter Lag (s)')
ylabel('Filter (norm.)')
title('Filter from log stimulus')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end



return


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

	