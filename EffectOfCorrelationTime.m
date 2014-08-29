%% Effect of Correlation Time on Linear Filters
% What is the effect of correlation time of the stimulus on the ORN response, and the filter that we calculate from this data? We have data where the same class of ORN responds to flickering stimuli with different correlation lengths: 30ms, 50ms and 100ms. The filters extracted from these datasets are slightly different, and fitted output non-linearities seem to have slopes that decrease with increasing correlation time. 

%%
% Is this a real phenomenon? Why does this happen?

load('/local-data/DA-paper/data.mat')

%%
% First, we plot the autocorrelation functions of the valve, the stimulus and the ORN and check that they make sense. 

dt = 0.003;
at = dt:dt:100*dt;
a30 = autocorr(data(4).Valve,99);
a100 = autocorr(data(3).Valve,99);
a50 = autocorr(data(5).Valve,99);

figure('outerposition',[0 0 1300 500],'PaperUnits','points','PaperSize',[1300 500]); hold on
subplot(1,3,1), hold on
plot(at,a30,'k')
plot(at,a50,'r')
plot(at,a100,'g')
set(gca,'XLim',[0 0.3])
xlabel('lag (s)')
ylabel('Autocorrelation')
title('Valve')
PrettyFig;
legend 30ms 50ms 100ms

a30 = autocorr(data(4).PID,99);
a100 = autocorr(data(3).PID,99);
a50 = autocorr(data(5).PID,99);

subplot(1,3,2), hold on
plot(at,a30,'k')
plot(at,a50,'r')
plot(at,a100,'g')
title('PID')
set(gca,'XLim',[0 0.3])
legend 30ms 50ms 100ms

a30 = autocorr(data(4).ORN,99);
a100 = autocorr(data(3).ORN,99);
a50 = autocorr(data(5).ORN,99);

subplot(1,3,3), hold on
plot(at,a30,'k')
plot(at,a50,'r')
plot(at,a100,'g')
set(gca,'XLim',[0 0.3])
title('ORN')
PrettyFig;
legend 30ms 50ms 100ms

%% Are filters different for different correlations?
% Now we fit filters to these data sets. In the figure below, the filters on the left are calculated using a fixed regularisation (of 1), and the filters on the right, the regularisation is varied in each case to find the "best" filter, minimising errors, high frequency components, and having a gain as close to unity as possible. In either case, we see the filters are slightly different for the three data sets; the correlation time of the stimulus seems to affect the filter we extract. 

K30=FindBestFilter(data(4).PID, data(4).ORN, [],'min_cutoff =0;','regmax =1;','regmin =1;');
K50=FindBestFilter(data(3).PID, data(3).ORN, [],'min_cutoff =0;','regmax =1;','regmin =1;');
[K100,~,filtertime]=FindBestFilter(data(5).PID, data(5).ORN, [],'min_cutoff =0;','regmax =1;','regmin =1;');

filtertime = filtertime*dt;

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
plot(filtertime,K30/max(K30),'k')
plot(filtertime,K50/max(K50),'r')
plot(filtertime,K100/max(K100),'g')
legend 30ms 50ms 100ms
set(gca,'XLim',[min(filtertime) max(data(4).filtertime)],'YLim',[-0.5 1.2])
title('Fixed regularisation')
xlabel('Filter Lag')
PrettyFig;

subplot(1,2,2), hold on
plot(data(4).filtertime,data(4).K/max(data(4).K),'k')
plot(data(3).filtertime,data(3).K/max(data(3).K),'r')
plot(data(5).filtertime,data(5).K/max(data(5).K),'g')
legend 30ms 50ms 100ms
set(gca,'XLim',[min(data(4).filtertime) max(data(4).filtertime)],'YLim',[-0.5 1.2])
title('"Best" regularisation')
xlabel('Filter Lag')
PrettyFig;

%%
% These filters differ not only in their shape but also their amplitude, which is a simple scaling because the input statistics are different. 

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(data(4).filtertime,data(4).K,'k')
plot(data(3).filtertime,data(3).K,'r')
plot(data(5).filtertime,data(5).K,'g')
legend 30ms 50ms 100ms
set(gca,'XLim',[min(data(4).filtertime) max(data(4).filtertime)])
xlabel('Filter Lag (s)')
PrettyFig;

%%
% Why does the filter amplitude change so much? let's look at the stimulus and the response histograms:

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
[y,x] = hist(data(4).PID);
plot(x,y,'k')
[y,x] = hist(data(3).PID);
plot(x,y,'r')
[y,x] = hist(data(5).PID);
plot(x,y,'g')
legend 30ms 50ms 100ms
xlabel('Stimulus Amplitude (PID, V)')
PrettyFig;

subplot(1,2,2), hold on
[y,x] = hist(data(4).ORN);
plot(x,y,'k')
[y,x] = hist(data(3).ORN);
plot(x,y,'r')
[y,x] = hist(data(5).ORN);
plot(x,y,'g')
legend 30ms 50ms 100ms
xlabel('ORN Response (Hz)')
PrettyFig;

%% Cross-prediction using different filters
% To disentangle the effects of filter amplitude from possible changes in filter shape, we fit a static non-linearity to the output of each linear filter and use that to predict the output.

this_filter = data(4).K/max(data(4).K);
fp = convolve(data(4).time,data(4).PID,this_filter,data(4).filtertime);
[~, NLN(1,1), NLN(1,2), NLN(1,3)] = FitNonLinearity(fp,data(4).ORN,'hill');

this_filter = data(3).K/max(data(3).K);
fp = convolve(data(3).time,data(3).PID,this_filter,data(3).filtertime);
[~, NLN(2,1), NLN(2,2), NLN(2,3)] = FitNonLinearity(fp,data(3).ORN,'hill');

this_filter = data(5).K/max(data(5).K);
fp = convolve(data(5).time,data(5).PID,this_filter,data(5).filtertime);
[~, NLN(3,1), NLN(3,2), NLN(3,3)] = FitNonLinearity(fp,data(5).ORN,'hill');

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ss = max(fp)/200;
x = min(fp):ss:max(fp);
plot(x,hill(NLN(1,:),x),'k')
plot(x,hill(NLN(2,:),x),'r')
plot(x,hill(NLN(3,:),x),'g')
xlabel('Input to NonLinearity')
ylabel('Function Output (Hz)')
legend 30ms 50ms 100ms
PrettyFig;

%%
% No monotonic relationship is observed between correlation length and slope of fitted non-linearity. 


%%
% How good are these filters at estimating the data? In the figure we below, we show the two extreme cases (the 30ms and 100ms correlated data) and the predictions from the filters calculated from their own data sets, along with predictions from filters calculated using the others' data. The cross-prediction filter outputs are then passed through the corresponding static non-linearities, 


figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
subplot(2,1,1), hold on
plot(data(4).time,data(4).ORN,'k')

% self predict
this_filter = data(4).K/max(data(4).K);
fp = convolve(data(4).time,data(4).PID,this_filter,data(4).filtertime);
fp = FitNonLinearity(fp,data(4).ORN,'hill');
plot(data(4).time,fp,'g')

% cross predict
this_filter = data(5).K/max(data(5).K);
fp = convolve(data(4).time,data(4).PID,this_filter,data(5).filtertime);
fp = FitNonLinearity(fp,data(4).ORN,'hill');
plot(data(4).time,fp,'r')

set(gca,'XLim',[mean(data(4).time)-3  mean(data(4).time)+2])
legend 30msData SamePrediction DifferentPrediction
ylabel('Firing Rate (Hz)')
PrettyFig;

subplot(2,1,2), hold on
plot(data(5).time,data(5).ORN,'k')

% self predict
this_filter = data(5).K/max(data(5).K);
fp = convolve(data(5).time,data(5).PID,this_filter,data(5).filtertime);
fp = FitNonLinearity(fp,data(5).ORN,'hill');
plot(data(5).time,fp,'g')

% cross predict
this_filter = data(4).K/max(data(4).K);
fp = convolve(data(5).time,data(5).PID,this_filter,data(5).filtertime);
fp = FitNonLinearity(fp,data(5).ORN,'hill');
plot(data(5).time,fp,'r')

set(gca,'XLim',[mean(data(5).time)-3  mean(data(5).time)+2])
legend 100msData SamePrediction DifferentPrediction
ylabel('Firing Rate (Hz)')
PrettyFig;

%%
% How good are these cross-predictions? For the three data sets we have, we measure how well the filter from each data set predicts the data. In the matrix below, the element in ith row and jth column shows the r-square of linear prediction by the ith filter to jth data set, where the first data set has 30ms correlated stimulus, the second has 50ms correlated stimulus, and the 3rd has 100ms correlated stimulus. 

r = NaN(3);
this_data = [4 3 5];

for i = 1:3
	for j = 1:3
		% use ith filter to predict j data
		a = this_data(i);
		b= this_data(j);
		% this_filter = data(a).K/max(data(a).K);
		% this_filter = this_filter*max(data(b).K);
		this_filter = data(a).K;
		fp = mean(data(b).ORN) + convolve(data(a).time,data(b).PID,this_filter,data(a).filtertime);
		fp(fp<0)=0;

		% compute r-square
		r(i,j) = rsquare(fp,data(b).ORN);
	end
	clear j
end
clear i

disp(r)

%%
% In particular, the 30ms filter outperforms all other filters, even the "native" fitlers of the 50ms and 100ms data. 

%% Synthetic data: filter and non-linearity estimation
% What effect does using stimuli with various correlation lengths have on our estimation of the filter? Here we use a "fake" filter to represent the neuron, and pass stimuli with various correlation lengths to it, and then back out the filter from the synthetic data output. 

% make a fake filter
clear p t
t=1:202;
p.theta=0.2000;
p.k=100;
p2.theta=0.400;
p2.k=100;
fake_K = (GammaDist(t,p) -  0.7*GammaDist(t,p2));
fake_K = fake_K/max(fake_K);
fake_K = fake_K*max(data(4).K);
clear p

% make fake orn outputs
fakeORN1 = convolve(data(4).time,data(4).PID,fake_K,data(4).filtertime);
fakeORN2 = convolve(data(3).time,data(3).PID,fake_K,data(3).filtertime);
fakeORN3 = convolve(data(5).time,data(5).PID,fake_K,data(5).filtertime);
p = data(5).PID(1:length(data(5).PID)/2);
pt = data(5).time(1:2:end);
pt = interp1(pt,p,data(5).time);
pt(end) = 0;
fakeORN4 = convolve(data(5).time,pt,fake_K,data(5).filtertime);

% add noise
fakeORN1 = fakeORN1 + randn(length(fakeORN1),1);
fakeORN2 = fakeORN2 + randn(length(fakeORN2),1);
fakeORN3 = fakeORN3 + randn(length(fakeORN3),1);
fakeORN4 = fakeORN4 + randn(length(fakeORN4),1);

% back out filters
fake_K1=FindBestFilter(data(4).PID, fakeORN1, []);
fake_K2=FindBestFilter(data(3).PID, fakeORN2, []);
fake_K3=FindBestFilter(data(5).PID, fakeORN3, []);
fake_K4=FindBestFilter(pt, fakeORN4, []);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
filtertime=mean(diff(data(4).time)):mean(diff(data(4).time)):mean(diff(data(4).time))*length(fake_K1);
filtertime = filtertime - 0.1*max(filtertime);
plot(filtertime,fake_K1/max(fake_K1),'k')
plot(filtertime,fake_K2/max(fake_K2),'r')
plot(filtertime,fake_K3/max(fake_K3),'g')
plot(filtertime,fake_K4/max(fake_K4),'b')
legend 30ms 50ms 100ms 200ms 
title('Reconstructed Filters')
set(gca,'XLim',[min(filtertime) max(filtertime)])
xlabel('Filter Lag (s)')
PrettyFig;


%%
% We see that the estimated filter gets wider for longer correlation lengths, even though we know that the actual filter is constant. Thus, trying to estimate filters of neurons with longer correlation times can back out longer-than-actual timescales in the filter. This is because the response of the system to fast time scales isn't being probed at all, or there is something wrong with the way I regularise. 

% now back out nonlinearities 
% fake_K1 = fake_K1/max(fake_K1);
% fake_K2 = fake_K2/max(fake_K2);
% fake_K3 = fake_K3/max(fake_K3);
% fake_K4 = fake_K4/max(fake_K4);

% fp1 = convolve(data(4).time,data(4).PID,fake_K1,filtertime);
% fp2 = convolve(data(3).time,data(3).PID,fake_K2,filtertime);
% fp3 = convolve(data(5).time,data(5).PID,fake_K3,filtertime);
% fp4 = convolve(data(5).time,pt,fake_K4,filtertime);

% clear NLN
% [~, NLN(1,1), NLN(1,2), NLN(1,3)] = FitNonLinearity(fp1,fakeORN1,'hill');
% [~, NLN(2,1), NLN(2,2), NLN(2,3)] = FitNonLinearity(fp2,fakeORN2,'hill');
% [~, NLN(3,1), NLN(3,2), NLN(3,3)] = FitNonLinearity(fp3,fakeORN3,'hill');
% [~, NLN(4,1), NLN(4,2), NLN(4,3)] = FitNonLinearity(fp4,fakeORN4,'hill');

% subplot(1,2,2), hold on
% ss = max(fp1)/200;
% x = min(fp1):ss:max(fp1);
% plot(x,hill(NLN(1,:),x),'k')
% plot(x,hill(NLN(2,:),x),'r')
% plot(x,hill(NLN(3,:),x),'g')
% plot(x,hill(NLN(4,:),x),'b')
% xlabel('Input to NonLinearity')
% xlabel('Function Output (Hz)')
% legend 30ms 50ms 100ms 200ms
% PrettyFig;




%% Version Info
% The file that generated this document is called:
disp(mfilename)

%%
% and its md5 hash is:
Opt.Input = 'file';
disp(DataHash(strcat(mfilename,'.m'),Opt))
