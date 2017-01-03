


pHeader;

%% Kinetics of LFP and firing rate 
% In this document, I analyse the kinetics of firing rate and LFP w.r.t the stimulus using cross correlation functions. 

%%
% In the following figure, I how the lags of the LFP and firing rate vary with stimulus background from the gaussian fluctuation data. (a) Normalised cross correlation computed from stimulus to LFP for the lowest dose (purple) and for the highest dose (yellow). (b) Raw LFP traces at odour offset show that LFP kinetics slow down with increasing mean stimulus (brighter colours). (c) Normalised cross correlation computed from stimulus to firing rate for the lowest dose (purple) and for the highest dose (yellow). No major change is visible. (d) Lags of firing rate (black) and LFP (red) estimated from peaks of cross correlation functions, as a function of the mean stimulus in that trial. Error bars are standard error of the mean. Cross correlations were estimated from 20 seconds of data in each trial, chunked into 1-second long fragments, and then averaged. 


% get the filter from the Gaussian stimuli 
clearvars -except being_published 
MSGdata = consolidateData2(getPath(dataManager,'93ba5d68174e3df9f462a1fc48c581da'));
MSGdata = cleanMSGdata(MSGdata);
v2struct(MSGdata)

clear p
p.A = 30;
p.adap_tau = 500;
p.K_n = 2;
p.K_tau = 100;

% generate responses for every trial
LFP_model = NaN*PID;
for i = 1:size(PID,2)
	textbar(i,size(PID,2))
	LFP_model(:,i) = asNL_euler(PID(:,i),p);
end


LFP_lags = NaN*paradigm;
LFP_model_lags = NaN*paradigm;
LFP_max_corr = NaN*paradigm;
firing_lags = NaN*paradigm;
firing_max_corr = NaN*paradigm;
a = 35e3+1; z = 55e3;
chunk_size = 1e3;

LFP_xcorr = NaN(chunk_size*2 - 1,20,length(paradigm));
LFP_model_xcorr = NaN(chunk_size*2 - 1,20,length(paradigm));

fA_xcorr = NaN(chunk_size*2 - 1,20,length(paradigm));

S_a = NaN(chunk_size*2 - 1,20,length(paradigm));

% compute LFP and firing rate lags
for i = 1:length(paradigm)
	S = PID(a:z,i);
	X = -raw_LFP(a:z,i);
	X_pred = LFP_model(a:z,i);
	R = fA(a:z,i);

	% reshape into chunks
	S = reshape(S,chunk_size,length(S)/chunk_size);
	X = reshape(X,chunk_size,length(X)/chunk_size);
	X_pred = reshape(X_pred,chunk_size,length(X_pred)/chunk_size);
	R = reshape(R,chunk_size,length(R)/chunk_size);

	S = bsxfun(@minus, S, mean(S));
	X = bsxfun(@minus, X, mean(X));
	X_pred = bsxfun(@minus, X_pred, mean(X_pred));
	R = bsxfun(@minus, R, mean(R));

	X_lag = NaN(chunk_size*2-1,size(S,2));
	X_pred_lag = NaN(chunk_size*2-1,size(S,2));
	R_lag = NaN(chunk_size*2-1,size(S,2));
	S_acorr = NaN(chunk_size*2-1,size(S,2));
	for j = 1:size(S,2)
		X_lag(:,j) = xcorr(X(:,j)/std(X(:,j)),S(:,j)/std(S(:,j)));
		X_pred_lag(:,j) = xcorr(X_pred(:,j)/std(X_pred(:,j)),S(:,j)/std(S(:,j)));
		R_lag(:,j) = xcorr(R(:,j)/std(R(:,j)),S(:,j)/std(S(:,j)));
		S_acorr(:,j) = xcorr(S(:,j)/std(S(:,j)),S(:,j)/std(S(:,j)));
	end
	X_lag = X_lag/chunk_size;
	LFP_xcorr(:,:,i) = X_lag;
	X_lag = mean(X_lag,2);
	[LFP_max_corr(i),loc] = max(X_lag);
	LFP_lags(i) = loc - 1e3;

	X_pred_lag = X_pred_lag/chunk_size;
	LFP_model_xcorr(:,:,i) = X_pred_lag;
	X_pred_lag = mean(X_pred_lag,2);
	[~,loc] = max(X_pred_lag);
	LFP_model_lags(i) = loc - 1e3;

	R_lag = R_lag/chunk_size;
	fA_xcorr(:,:,i) = R_lag;
	R_lag = mean(R_lag,2);
	[firing_max_corr(i),loc] = max(R_lag);
	firing_lags(i) = loc - 1e3;

	S_acorr = S_acorr/chunk_size;
	S_a(:,:,i) = S_acorr;
end


firing_lags(firing_lags<-100) = NaN; % obviously wrong


figure('outerposition',[0 0 1501 500],'PaperUnits','points','PaperSize',[1501 500]); hold on
for i = 1:3
	ax(i) = subplot(1,3,i); hold on
end
ax(4) = ax(3);

set(ax(1),'XLim',[0 1],'YLim',[-.5 1.1])
set(ax(2),'XLim',[0 1],'YLim',[-.5 1.1])
xlabel(ax(1),'Lag (s)')
xlabel(ax(2),'Lag (s)')
ylabel(ax(1),'Cross correlation (norm)')
ylabel(ax(2),'Cross correlation (norm)')
title(ax(1),'Stimulus \rightarrow LFP')
title(ax(2),'Stimulus \rightarrow Firing rate')

set(ax(3),'XLim',[54 57])
xlabel(ax(3),'Time (s)')
ylabel(ax(3),'\DeltaLFP (norm)')

% compare LFP lags at the lowest and highest dose
c = parula(max(paradigm)+1);
temp = LFP_xcorr(:,:,paradigm==1);
temp = reshape(temp,2*chunk_size - 1, size(temp,2)*size(temp,3));
lags = 1e-3*(1:length(temp)) - 1;
plot(ax(1),lags,mean(temp,2)/max(mean(temp,2)),'Color',c(1,:))
temp = LFP_xcorr(:,:,paradigm==10);
temp = reshape(temp,2*chunk_size - 1, size(temp,2)*size(temp,3));
lags = 1e-3*(1:length(temp)) - 1;
plot(ax(1),lags,mean(temp,2)/max(mean(temp,2)),'Color',c(10,:))


% compare firing lags at the lowest and highest dose
c = parula(max(paradigm)+1);
temp = fA_xcorr(:,:,paradigm==1);
temp = reshape(temp,2*chunk_size - 1, size(temp,2)*size(temp,3));
lags = 1e-3*(1:length(temp)) - 1;
plot(ax(2),lags,mean(temp,2)/max(mean(temp,2)),'Color',c(1,:))
temp = fA_xcorr(:,:,paradigm==10);
temp = reshape(temp,2*chunk_size - 1, size(temp,2)*size(temp,3));
lags = 1e-3*(1:length(temp)) - 1;
plot(ax(2),lags,nanmean(temp,2)/max(nanmean(temp,2)),'Color',c(10,:))

% show lags of firing rate and LFP for all doses
mean_stim = mean(PID(a:z,:));
c = lines(10);
LFP_color = c(5,:);
firing_color = c(4,:);
for i = 1:max(paradigm)
	x = nonnans(firing_lags(paradigm==i));
	errorbar(ax(4),mean(mean_stim(paradigm==i)),mean(x),sem(x),'Color',LFP_color)
	x = nonnans(LFP_lags(paradigm==i));
	errorbar(ax(4),mean(mean_stim(paradigm==i)),mean(x),sem(x),'Color',firing_color)
end
xlabel(ax(4),'\mu_{Stimulus} (V)')
ylabel(ax(4),'Lag (ms)')
legend({'Firing rate','LFP'},'Location','northwest')
set(ax(4),'XLim',[0 1.7],'YLim',[0 200])

za = 54e3; zz = 57e3;
time = 1e-3*(1:length(PID));

for i = 1:max(paradigm)
	X = LFP(:,paradigm==i);
	for j = 1:width(X)
		X(:,j) = X(:,j) - mean(X(za:55e3,j));
	end
	X = nanmean(X,2);

	X = X/max(X(za:zz));
	plot(ax(3),time(za:zz),X(za:zz,:),'Color',c(i,:))

end

prettyFig;

labelFigure('column_first',true)

if being_published	
	snapnow	
	delete(gcf)
end

%% 
% also check that the autocorrelation functions of all the stimulus is the same
figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
c = parula(max(paradigm)+1);
for i = 1:max(paradigm)
	temp = S_a(:,:,paradigm==i);
	temp = reshape(temp,2*chunk_size - 1, size(temp,2)*size(temp,3));
	lags = 1e-3*(1:length(temp)) - 1;
	plot(lags,mean(temp,2)/max(mean(temp,2)),'Color',c(i,:))
end

xlabel('Lag (s)')
ylabel('Autocorrelation (norm)')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%% 
% Can we reproduce this slowdown with a model? 




%%
% Can we also see signs of this adaptation in kinetics in naturalistic odor stimuli? 

% % get this data
% clearvars -except being_published
% load('/local-data/DA-paper/data-for-paper/nat-stim/ab3A_nat_stim.ORNData','-mat')
% od(1) = [];
% for i = length(od):-1:1
% 	data(i).LFP = nanmean(od(i).LFP,2);
% 	data(i).stimulus = nanmean(od(i).stimulus,2);
% 	data(i).firing_rate = nanmean(od(i).firing_rate,2);
% end

% figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on

% X = data(1).firing_rate;

% for i = 1:length(loc)
% 	temp = X(loc(i)-200:loc(i)+500);
% 	temp = temp - mean(temp(1:200));
% 	plot(temp/max(temp(200:end)))
% end


% [~,loc] = findpeaks(S,'MinPeakProminence',.1);
% % ignore all peaks with another peak 200 ms afterwards
% rm_this = false(length(loc),1);
% for i = 1:length(loc)
% 	if any(find(loc<loc(i) + 300 & loc > loc(i)))
% 		rm_this(i) = true;
% 	end
% end

%% Version Info
%
pFooter;




