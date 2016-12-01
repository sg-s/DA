


pHeader;

%% Kinetics of LFP and firing rate 
% In this document, I analyse the kinetics of firing rate and LFP w.r.t the stimulus using cross correlation functions. 


% get the filter from the Gaussian stimuli 
clearvars -except being_published 
MSGdata = consolidateData2(getPath(dataManager,'93ba5d68174e3df9f462a1fc48c581da'));
MSGdata = cleanMSGdata(MSGdata);
v2struct(MSGdata)


LFP_lags = NaN*paradigm;
LFP_max_corr = NaN*paradigm;
firing_lags = NaN*paradigm;
firing_max_corr = NaN*paradigm;
a = 35e3+1; z = 55e3;
chunk_size = 1e3;

LFP_xcorr = NaN(chunk_size*2 - 1,20,length(paradigm));
fA_xcorr = NaN(chunk_size*2 - 1,20,length(paradigm));

% compute LFP and firing rate lags
for i = 1:length(paradigm)
	S = PID(a:z,i);
	X = -raw_LFP(a:z,i);
	R = fA(a:z,i);

	% reshape into chunks
	S = reshape(S,chunk_size,length(S)/chunk_size);
	X = reshape(X,chunk_size,length(X)/chunk_size);
	R = reshape(R,chunk_size,length(R)/chunk_size);

	S = bsxfun(@minus, S, mean(S));
	X = bsxfun(@minus, X, mean(X));
	R = bsxfun(@minus, R, mean(R));

	X_lag = NaN(chunk_size*2-1,size(S,2));
	R_lag = NaN(chunk_size*2-1,size(S,2));
	clear X_lag R_lag
	for j = 1:size(S,2)
		X_lag(:,j) = xcorr(X(:,j)/std(X(:,j)),S(:,j)/std(S(:,j)));
		R_lag(:,j) = xcorr(R(:,j)/std(R(:,j)),S(:,j)/std(S(:,j)));
	end
	X_lag = X_lag/chunk_size;
	LFP_xcorr(:,:,i) = X_lag;
	X_lag = mean(X_lag,2);
	[LFP_max_corr(i),loc] = max(X_lag);
	LFP_lags(i) = loc - 1e3;

	R_lag = R_lag/chunk_size;
	fA_xcorr(:,:,i) = R_lag;
	R_lag = mean(R_lag,2);
	[firing_max_corr(i),loc] = max(R_lag);
	firing_lags(i) = loc - 1e3;

end


firing_lags(firing_lags<-100) = NaN; % obviously wrong





figure('outerposition',[0 0 901 802],'PaperUnits','points','PaperSize',[901 802]); hold on
clear ax
for i = 1:4
	ax(i) = subplot(2,2,i); hold on
end
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
for i = 1:max(paradigm)
	x = nonnans(firing_lags(paradigm==i));
	errorbar(ax(4),mean(mean_stim(paradigm==i)),mean(x),sem(x),'k')
	x = nonnans(LFP_lags(paradigm==i));
	errorbar(ax(4),mean(mean_stim(paradigm==i)),mean(x),sem(x),'r')
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


if being_published	
	snapnow	
	delete(gcf)
end


%% Version Info
%
pFooter;




