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
% make a linear prediction using a filter fit to the mean data (this is almost exactly the same)
[K, ~, filtertime_full] = FindBestFilter(mean2(PID),mean2(fA),[],'regmax=1;','regmin=1;','filter_length=1999;','offset=500;');
filtertime_full = filtertime_full*mean(diff(tA));
filtertime = 1e-3*(-200:900);
K = interp1(filtertime_full,K,filtertime);
K = K/max(K);
fp = convolve(tA,mean2(PID),K,filtertime);
R = mean2(fA);
temp =fit(fp(~(isnan(fp) | isnan(R))),R(~(isnan(fp) | isnan(R))),'poly1');
fp = fp*temp.p1;
fp = fp+temp.p2;

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
plot(tA,mean2(fA)./fp,'r')

% cache the gain to use later
cache('gain',mean2(fA)./fp)
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
set(gca,'YLim',[.3 5],'YScale','log')

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
% If we exclude points where there is a poor fit to the data ($r^2 < .8$), and take the absolute value of the gain, we get the following picture:

gain(r2<.8) = NaN;
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


history_lengths = logspace(-1,1,30);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
clear ph
ph(3) = subplot(1,2,1);
ph(4) = subplot(1,2,2);

[p,~,~,~,~,history_lengths]=GainAnalysisWrapper2('response',mean2(fA),'prediction',fp,'stimulus',mean2(PID),'time',tA,'history_lengths',history_lengths,'example_history_length',history_lengths(11),'use_cache',1,'engine',@GainAnalysis5,'ph',ph);
set(ph(4),'XLim',[.09 11]) % to show .1 and 10 on the log scale


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
