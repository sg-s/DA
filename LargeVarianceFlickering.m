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

%          ######  ######## #### ##     ## ##     ## ##       ##     ##  ######  
%         ##    ##    ##     ##  ###   ### ##     ## ##       ##     ## ##    ## 
%         ##          ##     ##  #### #### ##     ## ##       ##     ## ##       
%          ######     ##     ##  ## ### ## ##     ## ##       ##     ##  ######  
%               ##    ##     ##  ##     ## ##     ## ##       ##     ##       ## 
%         ##    ##    ##     ##  ##     ## ##     ## ##       ##     ## ##    ## 
%          ######     ##    #### ##     ##  #######  ########  #######   ######  


figure('outerposition',[0 0 1400 1000],'PaperUnits','points','PaperSize',[1400 1000]); hold on
subplot(2,8,9:14), hold on
plot(tA,mean2(fA),'k')
set(gca,'XLim',[10 60])
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')

subplot(2,8,15:16), hold on
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

%        ########  ########  ######  ########   #######  ##    ##  ######  ######## 
%        ##     ## ##       ##    ## ##     ## ##     ## ###   ## ##    ## ##       
%        ##     ## ##       ##       ##     ## ##     ## ####  ## ##       ##       
%        ########  ######    ######  ########  ##     ## ## ## ##  ######  ######   
%        ##   ##   ##             ## ##        ##     ## ##  ####       ## ##       
%        ##    ##  ##       ##    ## ##        ##     ## ##   ### ##    ## ##       
%        ##     ## ########  ######  ##         #######  ##    ##  ######  ######## 


subplot(2,8,1:6), hold on
plot(tA,mean2(PID),'k')
set(gca,'XLim',[10 60])
xlabel('Time (s)')
ylabel('Odor Concentration (V)')

subplot(2,8,7:8), hold on
hash = DataHash(PID);
cached_data = cache(hash);
if isempty(cached_data)
	r2 = rsquare(PID);
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


figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1), hold on
for i = 1:width(PID)
	[y,x] = histcounts(PID(:,i),300);x(1) = [];
	plot(x,y)
end

set(gca,'YScale','log')
xlabel('PID (V)')
ylabel('Count')
title('Stimulus Histogram')

subplot(1,3,2), hold on
a = [];
for i = 1:width(PID)
	[zc,temp,c]=FindCorrelationTime(PID(:,i));
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
	[zc,temp,c]=FindCorrelationTime(fA(:,i));
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

%     ##       #### ##    ## ########    ###    ########     ##     ##  #######  ########  ######## ##       
%     ##        ##  ###   ## ##         ## ##   ##     ##    ###   ### ##     ## ##     ## ##       ##       
%     ##        ##  ####  ## ##        ##   ##  ##     ##    #### #### ##     ## ##     ## ##       ##       
%     ##        ##  ## ## ## ######   ##     ## ########     ## ### ## ##     ## ##     ## ######   ##       
%     ##        ##  ##  #### ##       ######### ##   ##      ##     ## ##     ## ##     ## ##       ##       
%     ##        ##  ##   ### ##       ##     ## ##    ##     ##     ## ##     ## ##     ## ##       ##       
%     ######## #### ##    ## ######## ##     ## ##     ##    ##     ##  #######  ########  ######## ######## 


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

%%
% What happens when we explicitly specify a filter that we obtained from stimulation with pseudo-Gaussian stimuli? How well does it perform? 

clear p
p.  tau1= 3*5.0122;
p.   K_n= 4.4572;
p.  tau2= 3*20.3750;
p.     A= 41.3682;
p.     n= 4.4297;
p.    Kd= 20.0591;
p.offset= 20.0278;
p.   K_A= 0.3491;

% save for later
Linear_Filter_p = p;

t = 1:1e3;
K = filter_gamma2(p.tau1,p.K_n,p.tau2,p.K_A,t);
t = t*mean(diff(tA));
fp = filter(K,1,mean2(PID));

% trivial scaling
fp = fp -0.6562;
fp = fp*1.7987;

fp_gaussian_K = fp;

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,4,1:3), hold on
plot(tA,mean2(fA),'k')
l=plot(tA,fp,'r');
r2 = rsquare(fp,mean2(fA));
legend(l,strcat('r^2=',oval(r2,3)))
xlabel('Time (s)')
ylabel('Firing rate (Hz)')

subplot(1,4,4), hold on
plot(t,K,'r')
xlabel('Filter Lag (s)')
ylabel('Filter (norm.)')
title('Filter from Gaussian Noise')


PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%  ######      ###    #### ##    ##       ###    ##    ##    ###    ##       ##    ##  ######  ####  ######  
% ##    ##    ## ##    ##  ###   ##      ## ##   ###   ##   ## ##   ##        ##  ##  ##    ##  ##  ##    ## 
% ##         ##   ##   ##  ####  ##     ##   ##  ####  ##  ##   ##  ##         ####   ##        ##  ##       
% ##   #### ##     ##  ##  ## ## ##    ##     ## ## ## ## ##     ## ##          ##     ######   ##   ######  
% ##    ##  #########  ##  ##  ####    ######### ##  #### ######### ##          ##          ##  ##        ## 
% ##    ##  ##     ##  ##  ##   ###    ##     ## ##   ### ##     ## ##          ##    ##    ##  ##  ##    ## 
%  ######   ##     ## #### ##    ##    ##     ## ##    ## ##     ## ########    ##     ######  ####  ######  


%% Gain Changes
% In this section, we look at gain changes in the ORN. We do so, first, by dividing the instantaneous firing rate by the linear prediction. 

figure('outerposition',[0 0 1000 450],'PaperUnits','points','PaperSize',[1000 450]); hold on
plot([-1 61],[1 1],'k--')
plot(tA,mean2(fA)./fp_normal,'r')
% a = mean2(fA);
% a(end) = [];
% b = fp;
% b(end) = [];
% a = mean(reshape(a,10,length(a)/10));
% b = mean(reshape(b,10,length(b)/10));
% plot(diff(a)./diff(b))

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
% It looks like the gain is changing around two fold in some points. In the following figure, we do a more detailed gain analysis, splitting the data according to when the stimulus is high or low in the past and checking the gain in those points (as before). 

% do gain analysis
clear x
x.response = mean2(fA); 
x.prediction = fp_normal;
x.stimulus = mean2(PID); 
x.time = tA;
x.filter_length = 299;
ph = [];

rm_this = [find(isnan(mean2(fA))) find(isnan(fp_normal)) ];
x.response(rm_this) = [];
x.prediction(rm_this) = [];
x.stimulus(rm_this) = [];
x.time(rm_this) = [];

history_lengths = (3*floor(1000*logspace(-1.5,1,30)/3))/1e3;
example_history_length = 0.135;

f2=figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ph(3) = subplot(1,2,1); hold on 
axis square
ph(4) = subplot(1,2,2); hold on

hash = DataHash(x);
cached_data = cache(hash);
if isempty(cached_data)
	[p_K,l,h] = GainAnalysis4(x,history_lengths,example_history_length,ph);
	cache(hash,p_K);
	% also cache the example history length
	s=abs(l-h);
	s(p_K(1,:)>0.05)=NaN;
	[~,loc]=max(s);
	ehl = history_lengths(loc);
	cache(DataHash(p_K),ehl);

else
	p_K = cached_data;
	ehl = cache(DataHash(p_K));
	GainAnalysis4(x,history_lengths,ehl,ph,p_K);
end

xlabel(ph(3),'Linear Prediction (Hz)')
set(ph(4),'XScale','log')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%%
% For completeness, we repeat the gain analysis, but this time for the prediction from the Gaussian Noise. 

% do gain analysis
clear x
x.response = mean2(fA); 
x.prediction = fp_gaussian_K;
x.stimulus = mean2(PID); 
x.time = tA;
x.filter_length = 299;
ph = [];

rm_this = [find(isnan(mean2(fA))) find(isnan(fp_gaussian_K)) ];
x.response(rm_this) = [];
x.prediction(rm_this) = [];
x.stimulus(rm_this) = [];
x.time(rm_this) = [];

history_lengths = (3*floor(1000*logspace(-1.5,1,30)/3))/1e3;
example_history_length = 0.135;

f2=figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ph(3) = subplot(1,2,1); hold on 
axis square
ph(4) = subplot(1,2,2); hold on

hash = DataHash(x);
cached_data = cache(hash);
if isempty(cached_data)
	[p_gK,l,h] = GainAnalysis4(x,history_lengths,example_history_length,ph);
	cache(hash,p_gK);
	% also cache the example history length
	s=abs(l-h);
	s(p_gK(1,:)>0.05)=NaN;
	[~,loc]=max(s);
	ehl = history_lengths(loc);
	cache(DataHash(p_gK),ehl);

else
	p_gK = cached_data;
	ehl = cache(DataHash(p_gK));
	GainAnalysis4(x,history_lengths,ehl,ph,p_gK);
end

xlabel(ph(3),'Linear Prediction (Hz)')
set(ph(4),'XScale','log')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end


%% Does an output nonlinearity explain all observed gain changes?  
% In this section, we fit an output nonlinearity to the linear prediction *post-hoc*, and repeat the gain analysis to see if this operation can account for all previously observed gain changes. The output nonlinearity is chosen to be a smoothing spline. 

y = mean2(fA); y = y(:);
cf =fit(fp(1:10:end),y(1:10:end),'smoothingspline','SmoothingParam',0.01); 
figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(fp,y,'.','Color',[.8 .8 .8]), hold on, plot(sort(fp),cf(sort(fp)))
xlabel('Filter Output (Hz)')
ylabel('Neuron Response (Hz)')
% cache for later use
cache('smoothingspline_LN',cf);

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end



% clear p
% p.A= 56.9845;
% p.k= 23.5616;
% p.n= 2.9577;

% fp_hill = hill(p,fp_normal);
fp_ss = cf(fp_normal);

% do gain analysis
clear x
x.response = mean2(fA); 
x.prediction = fp_ss;
x.stimulus = mean2(PID); 
x.time = tA;
x.filter_length = 299;
ph = [];

rm_this = [find(isnan(mean2(fA))) find(isnan(fp_ss)) ];
x.response(rm_this) = [];
x.prediction(rm_this) = [];
x.stimulus(rm_this) = [];
x.time(rm_this) = [];

history_lengths = (3*floor(1000*logspace(-1.5,1,30)/3))/1e3;
example_history_length = 0.135;

f2=figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ph(3) = subplot(1,2,1); hold on 
axis square
ph(4) = subplot(1,2,2); hold on

hash = DataHash(x);
cached_data = cache(hash);
if isempty(cached_data)
	[p_LN,l,h] = GainAnalysis4(x,history_lengths,example_history_length,ph);
	cache(hash,p_LN);
	% also cache the example history length
	s=abs(l-h);
	s(p_LN(1,:)>0.05)=NaN;
	[~,loc]=max(s);
	ehl = history_lengths(loc);
	cache(DataHash(p_LN),ehl);

else
	p_LN = cached_data;
	ehl = cache(DataHash(p_LN));
	GainAnalysis4(x,history_lengths,ehl,ph,p_LN);
end

xlabel(ph(3),'LN Prediction (Hz)')
set(ph(4),'XScale','log')

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

