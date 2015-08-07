% Natural Flickering.m
% analysis of responses to very widely distributed stimuli
% 
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


%% Natural Flickering
% In this document, we generate odor stimuli flickers over a large range, to mimic what the fly might encounter in the "real" world. The idea is that odour stimuli are very broadly distributed in intensity, and that odour stimuli arrive in whiffs and clumps of whiffs (like in Vergassola et al.). 

%% Raw Data
% In the following data we look at the raw data from this experiment. On the left are time series of the stimulus (top, ethyl acetate) and the response of the ab3A neuron. On the right, we show the pair-wise coefficient of determination between all trials, and the pairwise slope of residuals between all trials. Values close to zero in both plots indicate that a) subsequent trials are well correlated and b) there is little trial-to-trial drift in the stimulus/response. 


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
% The following figure shows the odor stimulus (ethyl acetate). From left to right, we see the stimulus histogram showing the very long-tailed distribution of odour stimuli, the autocorrelation function of the stimulus, and the finally the autocorrelation function of the response. 


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
% In this section we fit a linear kernel from the stimulus to the response. 

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[900 500]); hold on
plot([-.1 1],[0 0 ],'k--')
[K, filtertime_full] = fitFilter2Data(mean2(PID),mean2(fA),'reg',1,'filter_length',1999,'offset',500);
filtertime_full = filtertime_full*mean(diff(tA));
filtertime = 1e-3*(-200:900);
K = interp1(filtertime_full,K,filtertime);

plot(filtertime,K,'r')
ylabel('Filter Amplitude')
xlabel('Filter Lag (s)')
title('Filter from stimulus')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end


%%
% We will use this filter to make a prediction of the response:

% convolve with filter to make prediction
fp = convolve(tA,mean2(PID),K,filtertime);

% correct for trivial scaling
R = mean2(fA);
temp =fit(fp(~(isnan(fp) | isnan(R))),R(~(isnan(fp) | isnan(R))),'poly1');
fp = fp*temp.p1;
fp = fp+temp.p2;


figure('outerposition',[0 0 1500 400],'PaperUnits','points','PaperSize',[1500 400]); hold on
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
title('Filter from stimulus')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end


%             #######  ##     ## ######## ########  ##     ## ######## 
%            ##     ## ##     ##    ##    ##     ## ##     ##    ##    
%            ##     ## ##     ##    ##    ##     ## ##     ##    ##    
%            ##     ## ##     ##    ##    ########  ##     ##    ##    
%            ##     ## ##     ##    ##    ##        ##     ##    ##    
%            ##     ## ##     ##    ##    ##        ##     ##    ##    
%             #######   #######     ##    ##         #######     ##    
           
%               ###    ##    ##    ###    ##       ##    ##  ######  ####  ######  
%              ## ##   ###   ##   ## ##   ##        ##  ##  ##    ##  ##  ##    ## 
%             ##   ##  ####  ##  ##   ##  ##         ####   ##        ##  ##       
%            ##     ## ## ## ## ##     ## ##          ##     ######   ##   ######  
%            ######### ##  #### ######### ##          ##          ##  ##        ## 
%            ##     ## ##   ### ##     ## ##          ##    ##    ##  ##  ##    ## 
%            ##     ## ##    ## ##     ## ########    ##     ######  ####  ######  

%% Output Analysis
% In this section, we analyse the prediction of the linear kernel in some more detail. The following figure shows a plot of the linear predictions vs. the actual response. We colour the data by the mean stimulus in the preceding 500ms.  We see that each excursion has a different slope, and it looks like the slopes correspond to the stimulus in the preceding 500ms. Larger mean stimulus in the preceding 500ms is indicated by brighter colours. In the second subplot, we plot the gain of each of these excursions vs the mean stimulus in the preceding 500ms. The red line is a power-law fit, and the exponent is mentioned in the fit legend. 

shat = ComputeSmoothedStimulus(mean2(PID),500);
% shat(fp<10) = 0;
shat = shat-min(shat);
shat = shat/max(shat);
shat = 1+ceil(shat*99);
shat(isnan(shat)) = 1;

ss = 1;
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on

cc = parula(100);
c= cc(shat,:);
scatter(fp(1:ss:end),R(1:ss:end),[],c(1:ss:end,:),'filled')

xlabel('Linear Prediction (Hz)')
ylabel('Actual response (Hz)')

subplot(1,2,2), hold on
fp = convolve(tA,mean2(PID),K,filtertime);
shat = ComputeSmoothedStimulus(mean2(PID),500);

% find all excursions (defined as firing rate crossing 10Hz)
[whiff_starts,whiff_ends] = ComputeOnsOffs(R>10);
mean_stim = NaN*whiff_ends;
gain = NaN*whiff_ends;
gain_err =  NaN*whiff_ends;
for i = 1:length(whiff_ends)
	mean_stim(i) = mean(shat(whiff_starts(i):whiff_ends(i)));
	ff=fit(fp(whiff_starts(i):whiff_ends(i)),R(whiff_starts(i):whiff_ends(i)),'poly1');
	gain(i) = ff.p1;
	temp = confint(ff);
	gain_err(i) = diff(temp(:,1))/2;
end
rm_this = (abs(gain_err./gain)) > .5; % throw out points where the estimate of gain has a more than 50% error
gain(rm_this) = [];
gain_err(rm_this) = [];
mean_stim(rm_this) = [];

errorbar(mean_stim,gain,gain_err,'k+')

ff = fit(mean_stim(:),gain(:),'power1','Weights',1./gain_err);

l = plot(sort(mean_stim),ff(sort(mean_stim)),'r');
L = strcat('y=\alpha x^{\beta}, \beta=',oval(ff.b),'. r^2=',oval(rsquare(ff(mean_stim),gain)));
legend(l,L)
xlabel('Mean Stimulus in preceding 500ms (V)')
ylabel('Gain (Hz/V)')

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


	