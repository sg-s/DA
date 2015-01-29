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

%%
% In this document, we generate odor stimuli flickers over a large range, like the "natural" stimuli, but never goes to zero, so that the neuron should never silence (allowing us to accurately follow its response). 

%%
% The odor stimulus and the response of the ORN looks like this:

load('/local-data/DA-paper/large-variance-flicker/2015_01_28_CS_ab3_2_EA.mat')
PID = data(4).PID;
time = 1e-4*(1:length(PID));
all_spikes = spikes(4).A;
load('/local-data/DA-paper/large-variance-flicker/2015_01_28_CS_ab3_3_EA.mat')
PID = vertcat(PID,data(4).PID);
all_spikes = vertcat(all_spikes,spikes(4).A);

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
% What does the stimulus distribution look like? 

[y,x] = histcounts(PID,300);x(1) = [];
figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(x,y,'k')
xlabel('PID (V)')
ylabel('Count')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end


%% Fitting A Linear Model to this Data
% Now we will attempt to find the best fit filter from the stimulus to the data. The following figure shows the filter backed out of this data, for different values of regularisation (scaled by the mean eigenvalue of the covariance matrix). 


figure('outerposition',[0 0 1400 500],'PaperUnits','points','PaperSize',[1400 500]); hold on
subplot(1,3,1), hold on
[K, ~, filtertime] = FindBestFilter(mean2(PID),mean2(fA),[],'regmax=0;','regmin=0;','filter_length=599;','offset=49;');
filtertime = filtertime*mean(diff(tA));
plot(filtertime,K)
title('reg = 0')
xlabel('Filter Lag (s)')


subplot(1,3,2), hold on
[K, ~, filtertime] = FindBestFilter(mean2(PID),mean2(fA),[],'regmax=0.1;','regmin=0.1;','filter_length=599;','offset=49;');
filtertime = filtertime*mean(diff(tA));
plot(filtertime,K)
title('reg = 0.1')
xlabel('Filter Lag (s)')

subplot(1,3,3), hold on
[K, ~, filtertime] = FindBestFilter(mean2(PID),mean2(fA),[],'regmax=1;','regmin=1;','filter_length=599;','offset=49;');
filtertime = filtertime*mean(diff(tA));
plot(filtertime,K)
title('reg = 1')
xlabel('Filter Lag (s)')

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end


%%
% There are some weird things here. First, the filter estimation doesn't seem to work very well (none of the filters look particularly clean.). Furthermore, the filter, when we can sort of pull it out, seems to have only 1 lobe, and no negative lobe. This is strange, as we expect that in the presence of background, the negative lobe of the filter would get larger. 

%%
% Even though this looks weird, we will use this filter to make a prediction of the response:

fp  =convolve(tA,mean2(PID),K,filtertime);
fp = fp + 47.76;
fp = fp*.4774;

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
% Now, we use the cross correlation function as a kernel and build a linear model around this. The following figure shows how well this linear model performs: 

K=(xc(tc>-0.1 & tc < 1));
filtertime = (tc(tc>-0.1 & tc < 1));
fp  =convolve(tA,mean2(PID),K,filtertime);

fp = fp +30.4745;
fp  =fp*0.7485;

fp_xcorr = fp;

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
title('Cross Correlation')

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%%
% In conclusion, we tried to get the filter in two different ways, and we see that the filter inexplicably lacks a negative lobe. 

%%
% What happens when we explicitly specify a filter with a negative lobe? For example, what happens if we use a filter extracted from Gaussian noise with this data? How well does it perform? 

p.  tau1= 3*5.0122;
p.   K_n= 4.4572;
p.  tau2= 3*20.3750;
p.     A= 41.3682;
p.     n= 4.4297;
p.    Kd= 20.0591;
p.offset= 20.0278;
p.   K_A= 0.3491;

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


%% Gain Changes
% In this section, we look at gain changes in the ORN. We do so, first, by dividing the instantaneous firing rate by the linear prediction. 

figure('outerposition',[0 0 1400 750],'PaperUnits','points','PaperSize',[1400 750]); hold on
subplot(2,1,1), hold on
plot([-1 61],[1 1],'k--')
plot(tA,mean2(fA)./fp_normal,'r')
ylabel('Gain')
title('from Rev. corr. filter')
set(gca,'YLim',[0 3])

subplot(2,1,2), hold on
plot([-1 61],[1 1],'k--')
plot(tA,mean2(fA)./fp_xcorr,'r')
ylabel('Gain')
title('from cross corr. filter')
xlabel('Time (s)')
set(gca,'YLim',[0 3])


PrettyFig;
if being_published
	snapnow
	delete(gcf)
end


return
 
%%
% It looks like the gain is changing around two fold in some points. In the following figure, we do a more detailed gain analysis, splitting the data according to when the stimulus is high or low in the past and checking the gain in those points (as before). 

% do gain analysis
clear x
x.response = fA; 
x.prediction = fp_xcorr;
x.stimulus = PID; 
x.time = tA;
x.filter_length = 299;
ph = [];

rm_this = [find(isnan(fA)) find(isnan(fp_normal)) ];
x.response(rm_this) = [];
x.prediction(rm_this) = [];
x.stimulus(rm_this) = [];

history_lengths = (3*floor(1000*logspace(-1.5,1,30)/3))/1e3;
example_history_length = 0.135;

f2=figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ph(3) = subplot(1,2,1); hold on 
axis square
ph(4) = subplot(1,2,2); hold on

if redo
	[p_LN,l,h] = GainAnalysis4(x,history_lengths,example_history_length,ph);
	s=abs(l-h);
	s(p_LN(1,:)>0.05)=NaN;
	[~,loc]=max(s);

	% save it for later
	ehl = history_lengths(loc);
	pb = p_LN;

else
	GainAnalysis4(x,history_lengths,ehl,ph,pb);
end

xlabel(ph(3),'Linear Prediction (Hz)')
set(ph(4),'XScale','log')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%% Analysis of MFC reliability
% A key issue seems to be that the stimulus is not very reproducibale. This may make the rest of the analysis rubbish. In the following figure, we generate this flickering stimulus with differnet switching times, and plot the reliability of the MFC. Each plot shows the r-square between all pairs of 10 trials for each switching time. 

%%
% In the following figure, the PID values for the feedback control on the Alicat MFC were P = 500, I = 0, D = 6000.

load('/local-data/DA-paper/MFC-calibration/Switching_Time_Reliability.mat')

figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
for i = 1:length(results)
	r= rsquare(results(i).data(2).MFC200);
	subplot(2,3,i), hold on
	imagescnan(r);
	caxis([0 1])
	colorbar
	axis image
	axis off
	title(strcat('switching time =',mat2str(1000*switching_time(i)),'ms'))
end

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

