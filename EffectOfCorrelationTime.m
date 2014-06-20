%% Effect of Correlation Time on Linear Filters
% What is the effect of correlation time of the stimulus on the ORN response, and the filter that we calculate from this data? 

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
PrettyFig;
legend 30ms 50ms 100ms

a30 = autocorr(data(4).PID,99);
a100 = autocorr(data(3).PID,99);
a50 = autocorr(data(5).PID,99);

subplot(1,3,2), hold on
plot(at,a30,'k')
plot(at,a50,'r')
plot(at,a100,'g')
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
PrettyFig;
legend 30ms 50ms 100ms

%%
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

set(gca,'XLim',[min(filtertime) max(data(4).filtertime)],'YLim',[-0.5 1.2])
xlabel('Filter Lag')
PrettyFig;

subplot(1,2,2), hold on
plot(data(4).filtertime,data(4).K/max(data(4).K),'k')
plot(data(3).filtertime,data(3).K/max(data(3).K),'r')
plot(data(5).filtertime,data(5).K/max(data(5).K),'g')

set(gca,'XLim',[min(data(4).filtertime) max(data(4).filtertime)],'YLim',[-0.5 1.2])
xlabel('Filter Lag')
PrettyFig;

%%
% How good are these filters at estimating the data? In the figure we below, we show the two extreme cases (the 30ms and 100ms correlated data) and the predictions from the filters calculated from their own data sets, along with predictions from filters calculated using the others' data. 


figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
subplot(2,1,1), hold on
plot(data(4).time,data(4).ORN,'k')
plot(data(4).time,data(4).LinearFit,'g')
% cross predict
fp = mean(data(4).ORN) + convolve(data(4).time,data(4).PID,data(5).K,data(5).filtertime);
fp(fp<0)=0;
plot(data(4).time,fp,'r')
set(gca,'XLim',[mean(data(4).time)-3  mean(data(4).time)+2])
legend 30msData SamePrediction DifferentPrediction
ylabel('Firing Rate (Hz)')
PrettyFig;

subplot(2,1,2), hold on
plot(data(5).time,data(5).ORN,'k')
plot(data(5).time,data(5).LinearFit,'g')
% cross predict
fp = mean(data(5).ORN) + convolve(data(5).time,data(5).PID,data(4).K,data(4).filtertime);
fp(fp<0)=0;
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
		fp = mean(data(b).ORN) + convolve(data(a).time,data(b).PID,data(a).K,data(a).filtertime);
		fp(fp<0)=0;

		% compute r-square
		r(i,j) = rsquare(fp,data(b).ORN);
	end
	clear j
end
clear i

disp(r)


%%
% In particular, predictions using the fastest filter outperform predictions even from filters from the same data set. 


