% LargeVarianceFlickering2.m
% 
% 
% created by Srinivas Gorur-Shandilya at 3:26 , 04 February 2015. Contact me at http://srinivas.gs/contact/
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
% From previous experiments, we see that the response of ORNs to largely varying stimuli that mimics the natural odor plumes is particularly interesting: in that no model we have can precisely account for the data, and in that, from other results, we expect to see a large variation of gain of the ORNs to these stimuli. This document is continued from Large Variance Flickering.m.  In this document we examine another dataset with supposedly shorter correlated stimuli.


load('/local-data/DA-paper/large-variance-flicker/2015_01_22_CS_F1_ab3_3_EtAc.mat')
PID = data(6).PID;
time = 1e-4*(1:length(PID));
all_spikes = spikes(6).A;

% A spikes --> firing rate
hash = DataHash(full(all_spikes));
cached_data = cache(hash);
if isempty(cached_data)
	fA = spiketimes2f(all_spikes,time);
	cache(hash,fA);
else
	fA = cached_data;
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
caxis([0 1])
colorbar
axis image
axis off
title(strcat('mean slope = ',oval(mean(s(~isnan(s))),2)))


%        ########  ########  ######  ########   #######  ##    ##  ######  ######## 
%        ##     ## ##       ##    ## ##     ## ##     ## ###   ## ##    ## ##       
%        ##     ## ##       ##       ##     ## ##     ## ####  ## ##       ##       
%        ########  ######    ######  ########  ##     ## ## ## ##  ######  ######   
%        ##   ##   ##             ## ##        ##     ## ##  ####       ## ##       
%        ##    ##  ##       ##    ## ##        ##     ## ##   ### ##    ## ##       
%        ##     ## ########  ######  ##         #######  ##    ##  ######  ######## 


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
caxis([0 1])
colorbar
axis image
axis off
title(strcat('mean slope = ',oval(mean(s(~isnan(s))),2)))

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

return

%% 
% Given the very poor reproducibility of this dataset, we will analyse this trial-by-trial. 

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
	[zc,~,c,temp]=FindCorrelationTime(PID(:,i));
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
	[zc,~,c,temp]=FindCorrelationTime(fA(:,i));
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
plot(1:length(r2),r2,'k+')
set(gca,'YLim',[0 1])
xlabel('Trial')
ylabel('r^2')
title('Prediction Fit Quality')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%% Output Analysis 
% In this section, we plot the prediction vs. the data in the space of the data. This serves two purposes. First, it makes apparent the shape of the output nonlinearity in the corresponding LN model. Second, it shows if the ORN's response trajectories contains loops, which is a tell-tale signature of dynamics that the filter does not account for. 

% account for some trivial scaling
clear p d
hash = DataHash(TrialFilters_fp);
cached_data = cache(hash);
if isempty(cached_data)
	for i = 1:width(fA)
		p(i).A = 0;
		p(i).B = 1;
		d.stimulus = TrialFilters_fp(:,i);
		d.response = fA(:,i);
		p(i) = FitModel2Data(@SubtractScale,d,p(i));
		TrialFilters_fp(:,i) = TrialFilters_fp(:,i) + p(i).A;
		TrialFilters_fp(:,i) = TrialFilters_fp(:,i)*p(i).B;
	end
	% cache
	cache(hash,p)
else
	p = cached_data;
	for i = 1:width(fA)
		TrialFilters_fp(:,i) = TrialFilters_fp(:,i) + p(i).A;
		TrialFilters_fp(:,i) = TrialFilters_fp(:,i)*p(i).B;
	end
end

% fit a nonlinearity 
clear p d
hash = DataHash(TrialFilters_fp);
cached_data = cache(hash);
if isempty(cached_data)
	for i = 1:width(fA)
		p(i).A = 100;
		p(i).k = 50;
		p(i).n = 2;
		d.stimulus = TrialFilters_fp(:,i);
		d.response = fA(:,i);
		p(i) = FitModel2Data(@hill,d,p(i));
		TrialFilters_fp(:,i) = hill(p(i),TrialFilters_fp(:,i));
	end
	% cache
	cache(hash,p)
else
	p = cached_data;
	for i = 1:width(fA)
		TrialFilters_LN(:,i) = hill(p(i),TrialFilters_fp(:,i));
	end
end


figure('outerposition',[0 0 1200 700],'PaperUnits','points','PaperSize',[1200 700]); hold on
a = [];
ss = 20; % sub sampling for plot
for i = 1:width(fA)
	a(i)=autoplot(width(fA),i); hold on
	plot(TrialFilters_LN(1:ss:end,i),fA(1:ss:end,i),'.','Color',[.6 .6 .6])
	title(strcat('Trial #',mat2str(i)))
	plot(sort(TrialFilters_fp(1:ss:end,i)),hill(p(i),sort(TrialFilters_fp(1:ss:end,i))),'r')
end 

PrettyFig('EqualiseX=1;','EqualiseY=1;');

if being_published
	snapnow
	delete(gcf)
end

%%
% How does the addition of this output nonlinearity improve the fit? The following plot shows the fit quality ($r^2$) for each trial 

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
clear l
r2 = NaN(width(fA),1);
for i = 1:width(TrialFilters)
	r2(i) = rsquare(TrialFilters_fp(:,i),fA(:,i));
end
plot(1:length(r2),r2,'k+');
r2 = NaN(width(fA),1);
for i = 1:width(TrialFilters)
	r2(i) = rsquare(TrialFilters_LN(:,i),fA(:,i));
end
plot(1:length(r2),r2,'r+')
legend('Linear Fit','LN Fit')
set(gca,'YLim',[0 1])
xlabel('Trial')
ylabel('r^2')
title('Prediction Fit Quality')

PrettyFig();

if being_published
	snapnow
	delete(gcf)
end


%% Gain Changes 
% In this section, we analyse the gain changes in each of these trials. For each trial, we perform the gain analysis as before, and each figure corresponds to one trial. Note that the fits to all the clouds of green points are very poor, invalidating this approach. 

for i = 1:width(PID)
	ph=GainAnalysisWrapper(fA(:,i),TrialFilters_LN(:,i),PID(:,i),tA);
	title(ph(4),strcat('Trial #',mat2str(i)))
end

%% 
% It is possible that there are fast gain changes, but our analysis techniques don't allow us to see what they are. (See <https://github.com/sg-s/DA/issues/112 this issue on Github>)

%% Dynamical Adaptation Model Fits (nothing to show)
% Perhaps the gain of the ORN is changing, even though we are unable to see it because the way we analyse it is poor. If this is the case, we should be able to fit a Dynamical Adaptation Model to this dataset and do a better job explaining the data. Is this the case? 


%%
% Unable to fit a DA Model to the data that is something other than a trivial filter. The model performs much worse than the rev. corr. filter. 





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