% LargeVarianceFlickering.m
% 
% created by Srinivas Gorur-Shandilya at 3:30 , 19 January 2015. Contact me at http://srinivas.gs/contact/
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

%% Response of ORNs to flickering odor stimuli with large variances
% From previous experiments, we see that the response of ORNs to largely varying stimuli that mimics the natural odor plumes is particularly interesting: in that no model we have can precisely account for the data, and in that, from other results, we expect to see a large variation of gain of the ORNs to these stimuli. 


load('/local-data/DA-paper/large-variance-flicker/2015_01_28_CS_ab3_2_EA.mat')
PID = data(4).PID;
time = 1e-4*(1:length(PID));
all_spikes = spikes(4).A;
B_spikes = spikes(4).B;
load('/local-data/DA-paper/large-variance-flicker/2015_01_28_CS_ab3_3_EA.mat')
PID = vertcat(PID,data(4).PID);
all_spikes = vertcat(all_spikes,spikes(4).A);
B_spikes = vertcat(B_spikes,spikes(4).B);

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


%        ########  ########  ######  ########   #######  ##    ##  ######  ######## 
%        ##     ## ##       ##    ## ##     ## ##     ## ###   ## ##    ## ##       
%        ##     ## ##       ##       ##     ## ##     ## ####  ## ##       ##       
%        ########  ######    ######  ########  ##     ## ## ## ##  ######  ######   
%        ##   ##   ##             ## ##        ##     ## ##  ####       ## ##       
%        ##    ##  ##       ##    ## ##        ##     ## ##   ### ##    ## ##       
%        ##     ## ########  ######  ##         #######  ##    ##  ######  ######## 




figure('outerposition',[0 0 1400 1000],'PaperUnits','points','PaperSize',[1400 1000]); hold on
subplot(2,9,10:14), hold on
plot(tA,mean2(fA),'k')
set(gca,'XLim',[10 60])
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


%          ######  ######## #### ##     ## ##     ## ##       ##     ##  ######  
%         ##    ##    ##     ##  ###   ### ##     ## ##       ##     ## ##    ## 
%         ##          ##     ##  #### #### ##     ## ##       ##     ## ##       
%          ######     ##     ##  ## ### ## ##     ## ##       ##     ##  ######  
%               ##    ##     ##  ##     ## ##     ## ##       ##     ##       ## 
%         ##    ##    ##     ##  ##     ## ##     ## ##       ##     ## ##    ## 
%          ######     ##    #### ##     ##  #######  ########  #######   ######  


subplot(2,9,1:5), hold on
plot(tA,mean2(PID),'k')
set(gca,'XLim',[10 60])
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
set(gca,'XLim',[10 60])
xlabel('Time (s)')
ylabel('Firing Rate (A) (Hz)')
set(gca,'YLim',[0 70])

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
set(gca,'XLim',[10 60])
xlabel('Time (s)')
ylabel('Firing Rate (B) (Hz)')
set(gca,'YLim',[0 70])

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


f=figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1), hold on
for i = 1:width(PID)
	[y,x] = histcounts(PID(:,i),300);x(1) = [];
	plot(x,y)
end

set(gca,'YScale','log')
xlabel('PID (V)')
ylabel('Count')
title('Stimulus Histogram')

figure(f);
plot_here= subplot(1,3,2); hold on
a = [];
for i = 1:width(PID)
	[zc,temp,c]=FindCorrelationTime(PID(:,i));
	a = [a temp];
	l=plot(plot_here,1:3e3,c(1:3e3));
end

xlabel(plot_here,'Lag (ms)')
ylabel(plot_here,'Autocorrelation')
legend(l,strcat('\tau=',oval(a),'ms'))
set(plot_here,'XScale','log')
title(plot_here,'Stimulus')

figure(f);
plot_here = subplot(1,3,3); hold on
a = [];
for i = 1:width(fA)
	[zc,temp,c]=FindCorrelationTime(fA(:,i));
	a = [a temp];
	l=plot(plot_here,1:3e3,c(1:3e3));
end

title(plot_here,'Response')
xlabel(plot_here,'Lag (ms)')
ylabel(plot_here,'Autocorrelation')
legend(l,strcat('\tau=',oval(a),'ms'))
set(plot_here,'XScale','log')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%      ##       #### ##    ## ########    ###    ########     ######## #### ######## 
%      ##        ##  ###   ## ##         ## ##   ##     ##    ##        ##     ##    
%      ##        ##  ####  ## ##        ##   ##  ##     ##    ##        ##     ##    
%      ##        ##  ## ## ## ######   ##     ## ########     ######    ##     ##    
%      ##        ##  ##  #### ##       ######### ##   ##      ##        ##     ##    
%      ##        ##  ##   ### ##       ##     ## ##    ##     ##        ##     ##    
%      ######## #### ##    ## ######## ##     ## ##     ##    ##       ####    ##    


%% Fitting A Linear Model to this Data
% Now we will attempt to find the best fit filter from the stimulus to the data. The following figure shows the "best" filter backed out of this data. 


figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[900 500]); hold on
plot([-.1 1],[0 0 ],'k--')
[K, ~, filtertime_full] = FindBestFilter(mean2(PID),mean2(fA),[],'regmax=1;','regmin=1;','filter_length=1999;','offset=500;');
filtertime_full = filtertime_full*mean(diff(tA));
filtertime = 1e-3*(-200:900);
K = interp1(filtertime_full,K,filtertime);

plot(filtertime,K,'r')
ylabel('Filter Amplitude')
xlabel('Filter Lag (s)')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end


%%
% We will use this filter to make a prediction of the response:

K = K/max(K);
fp  =convolve(tA,mean2(PID),K,filtertime);

fp = fp + 20.1314;
fp = fp*1.1323;

fp_normal = fp;

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,4,1:3), hold on
plot(tA,mean2(fA),'k')
l=plot(tA,fp,'r');
r2 = rsquare(fp,mean2(fA));
legend(l,strcat('r^2=',oval(r2,2)))
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')

subplot(1,4,4), hold on
plot(filtertime,K,'r')
xlabel('Filter Lag (s)')
ylabel('Filter (norm.)')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%%
% The following plot shows a neuron-wise analysis of the model fit a la Geffen and Meister 2009. Points above the diagonal correspond to model fits to neurons that are exceed the variability in the data. 

% make a geffen-meister plot
figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
[qx, qy] = GeffenMeister(fA(1e4:end,1:10),fp(1e4:end));
plot(qx,qy,'k+')
[qx, qy] = GeffenMeister(fA(1e4:end,11:end),fp(1e4:end));
plot(qx,qy,'k+')

plot([0 6],[0 6],'k--')
set(gca,'XLim',[0 6],'YLim',[0 6])
xlabel('(P_{S}/P_{N})^{1/2}','interpreter','tex')
ylabel('(P_{S}/P_{R})^{1/2}','interpreter','tex')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end



%%
% To figure out what's going on here, we will inspect the cross correlation function between the stimulus and the response: 

[xc,tc]=xcorr(mean2(fA)-mean(mean2(fA)),mean2(PID)-mean(mean2(PID)));
xc = xc/max(xc);
tc = tc*mean(diff(tA));

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,4,1:3), hold on
plot(tc,xc)
xlabel('Lag (s) (PID \rightarrow Firing rate)')
ylabel('Cross-correlation')

subplot(1,4,4), hold on
plot(tc,xc)
set(gca,'XLim',[-1 1])
xlabel('Lag (s)')

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end


%%
% The Firing rate lags the PID trace by 66ms. 


%  ######      ###    #### ##    ## 
% ##    ##    ## ##    ##  ###   ## 
% ##         ##   ##   ##  ####  ## 
% ##   #### ##     ##  ##  ## ## ## 
% ##    ##  #########  ##  ##  #### 
% ##    ##  ##     ##  ##  ##   ### 
%  ######   ##     ## #### ##    ## 
%    ###    ##    ##    ###    ##       ##    ##  ######  ####  ######  
%   ## ##   ###   ##   ## ##   ##        ##  ##  ##    ##  ##  ##    ## 
%  ##   ##  ####  ##  ##   ##  ##         ####   ##        ##  ##       
% ##     ## ## ## ## ##     ## ##          ##     ######   ##   ######  
% ######### ##  #### ######### ##          ##          ##  ##        ## 
% ##     ## ##   ### ##     ## ##          ##    ##    ##  ##  ##    ## 
% ##     ## ##    ## ##     ## ########    ##     ######  ####  ######  



%% Gain Changes
% In this section, we look at gain changes in the ORN. We do so, first, by dividing the instantaneous firing rate by the linear prediction. 

figure('outerposition',[0 0 1000 450],'PaperUnits','points','PaperSize',[1000 450]); hold on
plot([-1 61],[1 1],'k--')
plot(tA,mean2(fA)./fp_normal,'r')

% cache the gain to use later
cache('gain',mean2(fA)./fp_normal)
ylabel('Gain')
set(gca,'YLim',[0 2.1],'XLim',[5 60])

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%%
% Another way to think about the gain is to define the gain as 
%
% $$ g(t)=\frac{dr}{ds} $$
%

%%
% In the following figure, we estimate the gain as so, by fitting lines to intervals 400ms long, and sliding along 10ms. 

[gain,r2,t] = EstimateGain(fp,mean2(fA),400,10);
t = t*1e-3;

figure('outerposition',[0 0 1400 700],'PaperUnits','points','PaperSize',[1400 700]); hold on
subplot(4,1,1:3)
plot(t,gain,'k')
ylabel('Gain')

subplot(4,1,4)
plot(t,r2,'k')
xlabel('Time (s)')
ylabel('r^2')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%%
% If we exclude points where there is a poor fit to the data ($r^2 < .5$), and take the absolute value of the gain, we get the following picture:

gain(r2<.5) = NaN;
figure('outerposition',[0 0 1400 500],'PaperUnits','points','PaperSize',[1400 500]); hold on
plot([0 60],[1 1],'k--')
plot(t,abs(gain),'r')
ylabel('Gain')
xlabel('Time (s)')
set(gca,'YScale','log')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end


%%
% It looks like the gain is changing around two fold in some points. In the following figure, we do a more detailed gain analysis, splitting the data according to when the stimulus is high or low in the past and checking the gain in those points (as before). 

ph=GainAnalysisWrapper(mean2(fA),fp_normal,mean2(PID),tA);
xlabel(ph(3),'Linear Prediction (Hz)')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end


%% Does an output nonlinearity explain all observed gain changes?  
% In this section, we fit an output nonlinearity to the linear prediction *post-hoc*, and repeat the gain analysis to see if this operation can account for all previously observed gain changes. The output nonlinearity is chosen to be a Hill function. 

% y = mean2(fA); y = y(:);

clear p
p.A = 57.2707;
p.k = 23.7470;
p.n = 2.9373;

fp_hill = hill(p,fp_normal);


% cf =fit([fp(1:10:end)],[y(1:10:end)],'smoothingspline','SmoothingParam',0.01); 
% cf =fit( fp(1:10:end), y(1:10:end),'smoothingspline','SmoothingParam',1e-2); 

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(fp_normal,mean2(fA),'.','Color',[.8 .8 .8]), hold on, plot(sort(fp_normal),hill(p,sort(fp_normal)))
xlabel('Filter Output (Hz)')
ylabel('Neuron Response (Hz)')
% cache for later use

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end



% fp_ss = cf(fp_normal);

ph=GainAnalysisWrapper(mean2(fA),fp_hill,mean2(PID),tA);
xlabel(ph(3),'LN Prediction (Hz)');

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%% Further analysis of filter shape
% We expect that the negative lobe of the filter to be more prominent in this dataset, as the neuron should be more differentiating, as the signal consists of fluctuations on top of a background. But we don't see that. Why? Perhaps because of the non-gaussian statistics of the input stimulus, the filter shape is screwed up. To check this, we generate synthetic data using the bilobed linear filter, and then back out a linear filter from this.

%%
% In the following figures, we try to back out a filter (red) from synthetic data constructed using the black filter. 

clear p
p = Linear_Filter_p;

t = 1:1e3;
K = filter_gamma2(p.tau1,p.K_n,p.tau2,p.K_A,t);
t = t*mean(diff(tA));
fp = filter(K,1,mean2(PID));
fp = fp + 2*randn(length(fp),1);

figure('outerposition',[0 0 1400 500],'PaperUnits','points','PaperSize',[1400 500]); hold on

subplot(1,3,1), hold on
plot(t,K,'k')
% now back out a filter from this
[K_back, ~, filtertime] = FindBestFilter(mean2(PID),fp,[],'regmax=.01;','regmin=.01;','filter_length=999;');
plot(filtertime*1e-3,K_back/max(K_back),'r')
title('reg=.01')
xlabel('Filter Lag (s)')
ylabel('Filter Amplitude (norm)')


subplot(1,3,2), hold on
plot(t,K,'k')
% now back out a filter from this
[K_back, ~, filtertime] = FindBestFilter(mean2(PID),fp,[],'regmax=.1;','regmin=.1;','filter_length=999;');
plot(filtertime*1e-3,K_back/max(K_back),'r')
title('reg=.1')
xlabel('Filter Lag (s)')


subplot(1,3,3), hold on
plot(t,K,'k')
% now back out a filter from this
[K_back, ~, filtertime] = FindBestFilter(mean2(PID),fp,[],'regmax=1;','regmin=1;','filter_length=999;');
plot(filtertime*1e-3,K_back/max(K_back),'r')
title('reg=1')
xlabel('Filter Lag (s)')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%%
% It looks like my filter estimation routines can do a good job even with this non-Gaussian stimulus statistics. 

%% Neuron-by-Neuron analysis
% This data set has trials from multiple neurons. Here, we group by neuron ID, and repeat the analysis, to see if there any differences between the neurons, and if this causes a spurious gain change. 

fA1 = mean2(fA(:,1:10));
fA2 = mean2(fA(:,11:20));
PID1 = mean2(PID(:,1:10));
PID2 = mean2(PID(:,11:20));

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(tA,fA1,'m')
plot(tA,fA2,'g');
legend('ORN 1','ORN 2')
r2 = rsquare(fA1,fA2);
title(strcat('r^2=',oval(r2,3)))
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%%
% They look very similar. How do the filters backed out of the two neurons look? In the following figure, we back out filters from the stimuli and the responses corresponding to each neuron separately, and compare them below. We also compare the neuron-wise response predictions. 


figure('outerposition',[0 0 1300 700],'PaperUnits','points','PaperSize',[1300 700]); hold on
clear l
l(1)=subplot(2,4,1:3); hold on
title('ORN 1')
plot(tA,fA1,'k')
ylabel('Firing Rate (Hz)')

l(2)=subplot(2,4,5:7); hold on
title('ORN 2')
plot(tA,fA2,'k')
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')


l(3)=subplot(4,4,[8 12]); hold on
plot([-.1 1],[0 0 ],'k--')

[K, ~, filtertime_full] = FindBestFilter(PID1,fA1,[],'regmax=1;','regmin=1;','filter_length=1999;','offset=500;');
filtertime_full = filtertime_full*mean(diff(tA));
filtertime = 1e-3*(-200:900);
K1 = interp1(filtertime_full,K,filtertime);
K1 = K1/max(K1);

plot(l(3),filtertime,K1,'m')

[K, ~, filtertime_full] = FindBestFilter(PID2,fA2,[],'regmax=1;','regmin=1;','filter_length=1999;','offset=500;');
filtertime_full = filtertime_full*mean(diff(tA));
filtertime = 1e-3*(-200:900);
K2 = interp1(filtertime_full,K,filtertime);
K2 = K2/max(K2);

plot(l(3),filtertime,K2,'g')
ylabel('Filter Amplitude')
xlabel('Filter Lag (s)')
title(strcat('r^2=',oval(rsquare(K1,K2))))

% make the predictions. 
fp1 = convolve(tA,PID1,K1,filtertime); 
fp2 = convolve(tA,PID2,K2,filtertime);

% trivial scaling
fp1 = fp1+19.5943;
fp1 = fp1*1.0654;

fp2 = fp2+20.6811;
fp2 = fp2*1.1950;

clear ll
ll(1) = plot(l(1),tA,fp1,'m');
ll(2) = plot(l(2),tA,fp2,'g');

r21 = rsquare(fp1,fA1);
r22 = rsquare(fp2,fA2);


legend(ll(1),strcat('r^2=',oval(r21,3)))
legend(ll(2),strcat('r^2=',oval(r22,3)))

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end


%%
% In the following section, we analyse the static output non-linearities of the neurons individually. 


% cf1 =fit([fp1(1:10:end-300)],[fA1(1:10:end-300)],'smoothingspline','SmoothingParam',0.01);
% cf2 =fit([fp2(1:10:end-300)],[fA2(1:10:end-300)],'smoothingspline','SmoothingParam',0.01);  

clear p 
p(1).A = 52.7523;
p(1).k = 21.6964;
p(1).n = 2.9652;

p(2).A = 61.2457;
p(2).k = 25.5119;
p(2).n = 2.9746;

figure('outerposition',[0 0 1300 700],'PaperUnits','points','PaperSize',[1300 700]); hold on
clear l
l(1)=subplot(2,4,1:3); hold on
title('ORN 1')
plot(tA,fA1,'k')
ylabel('Firing Rate (Hz)')

l(2)=subplot(2,4,5:7); hold on
title('ORN 2')
plot(tA,fA2,'k')
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')


l(3)=subplot(4,4,[8 12]); hold on
plot(1:90,hill(p(1),1:90),'m')
plot(1:90,hill(p(2),1:90),'g')

clear ll
ll(1)=plot(l(1),tA,hill(p(1),fp1),'m');
ll(2)=plot(l(2),tA,hill(p(2),fp2),'g');

r21 = rsquare(hill(p(1),fp1),fA1);
r22 = rsquare(hill(p(2),fp2),fA2);

legend(ll(1),strcat('r^2=',oval(r21,3)))
legend(ll(2),strcat('r^2=',oval(r22,3)))

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%%
% Now for the crucial bit. Do we still observe gain changes when we analyse the ORNs one by one? 

ph=GainAnalysisWrapper(fA1,hill(p(1),fp1),PID1,tA);
ylabel(ph(3),'ORN 2 Firing Rate (Hz)')
PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

% now do ORN 2
ph=GainAnalysisWrapper(fA2,hill(p(2),fp2),PID2,tA);
ylabel(ph(3),'ORN 2 Firing Rate (Hz)')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%% Trial by Trial Analysis
% In this section, we perform a detailed analysis of the data on a trial-by-trial basis. The reasons are two-fold. First, it is to check if the gain changes we see are apparent even in these trial-to-trial analyses. Second, it is to see if the data lends itself to analysis on an individual trial basis. If this is the case, then the range of potential experiments we can do dramatically increases, as does the sort of data we can analyse. 

TrialFilters = NaN(1101,width(fA));
TrialFilters_fp = fA;
for i = 1:width(fA)
	[K, ~, filtertime_full] = FindBestFilter(PID(:,i),fA(:,i),[],'regmax=1;','regmin=1;','filter_length=1999;','offset=500;');
	filtertime_full = filtertime_full*mean(diff(tA));
	filtertime = 1e-3*(-200:900);
	K = interp1(filtertime_full,K,filtertime);
	TrialFilters(:,i) = K/max(K);
end

figure('outerposition',[0 0 1400 500],'PaperUnits','points','PaperSize',[1400 500]); hold on
subplot(1,3,1), hold on
plot(filtertime,TrialFilters)
xlabel('Lag (s)')
ylabel('Filter Amplitude')

subplot(1,3,2), hold on
r2 = rsquare(TrialFilters);
imagescnan(r2)
caxis([0 1])
axis image
colorbar
axis off
title(strcat('Filter mean r^2 = ',oval(mean(r2(~isnan(r2))),2)))


subplot(1,3,3), hold on
% make all the predictions
r2 = NaN(width(fA),1);
for i = 1:width(TrialFilters)
	TrialFilters_fp(:,i) = convolve(tA,PID(:,i),TrialFilters(:,i),filtertime);
	r2(i) = rsquare(TrialFilters_fp(:,i),fA(:,i));
end
plot(1:20,r2,'k+')
set(gca,'YLim',[0 1])
xlabel('Trial')
ylabel('r^2')
title('Prediction Fit Quality')

PrettyFig;

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
disp(DataHash(strcat(mfilename,'.m'),Opt))

%%
% This file should be in this commit:
[status,m]=unix('git rev-parse HEAD');
if ~status
	disp(m)
end
