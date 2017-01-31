


pHeader;

%% Kinetics of LFP and firing rate 
% In this document, I analyse the kinetics of firing rate and LFP w.r.t the stimulus using cross correlation functions. 

%%
% In the following figure, I how the lags of the LFP and firing rate vary with stimulus background from the gaussian fluctuation data. (a) Normalised cross correlation computed from stimulus to LFP for the lowest dose (purple) and for the highest dose (yellow). (b) Raw LFP traces at odour offset show that LFP kinetics slow down with increasing mean stimulus (brighter colours). (c) Normalised cross correlation computed from stimulus to firing rate for the lowest dose (purple) and for the highest dose (yellow). No major change is visible. (d) Lags of firing rate (black) and LFP (red) estimated from peaks of cross correlation functions, as a function of the mean stimulus in that trial. Error bars are standard error of the mean. Cross correlations were estimated from 20 seconds of data in each trial, chunked into 1-second long fragments, and then averaged. 

c = lines(5);
LFP_color = c(4,:);
firing_color = c(5,:);


clear p
p.      E0 = 20.5000;
p.      E1 = 2.2500;
p. w_minus = 5.0625;
p.  w_plus = 6.5000;
p.tau_adap = 7.0625;
p.      kT = 35.9375;
p.   K_tau = 29.7500;




 p(2).      E0 = 56.5000;
 p(2).      E1 = 2.3750;
 p(2). w_minus = 6.0625;
 p(2).  w_plus = 10;
 p(2).tau_adap = 2.0625;
 p(2).      kT = 35.9375;
 p(2).   K_tau = 67.7500;


p(4).      E0 = 36.5000;
p(4).      E1 = 1.7500;
p(4). w_minus = 10.0625;
p(4).  w_plus = 6.5000;
p(4).tau_adap = 2.6875;
p(4).      kT = 75.9375;
p(4).   K_tau = 45.7500;


data_hashes = {'93ba5d68174e3df9f462a1fc48c581da','bcd4cf4fe12817d084a2b06f981161ee','cd6753c0e4cf02895cd5e2c5cb58aa1a','3ea08ccfa892c6545d74bbdaaa6cbee1'};
odour_names = {'ethyl-acetate','1-pentanol','1-pentanol','2-butanone'};
orn_names = {'ab3A','ab3A','ab2A','ab2A'};
x_limits = [1.7 .03 .16 1.5];

for di = 1:length(data_hashes)

	MSGdata = consolidateData2(getPath(dataManager,data_hashes{di}));
	MSGdata = cleanMSGdata(MSGdata,'extract_filter',false);
	v2struct(MSGdata)

	% remove baseline from PID for each trial
	PID = bsxfun(@minus,PID,min(PID));


	% generate responses for every trial
	try
		LFP_model = NaN*PID;
		for i = 1:size(PID,2)
			LFP_model(30e3:55e3,i) = asNL3(PID(30e3:55e3,i),p(di));
			%LFP_model(:,i) = asNL_euler(PID(:,i),p(di));
		end
	catch
	end


	LFP_lags = NaN*paradigm;
	LFP_model_lags = NaN*paradigm;
	S_autocorr = NaN*paradigm;
	LFP_max_corr = NaN*paradigm;
	firing_lags = NaN*paradigm;
	firing_max_corr = NaN*paradigm;
	a = 30e3+1; z = 50e3;
	chunk_size = 1e3;

	LFP_xcorr = NaN(chunk_size*2 - 1,20,length(paradigm));

	fA_xcorr = NaN(chunk_size*2 - 1,20,length(paradigm));

	S_a = NaN(chunk_size*2 - 1,20,length(paradigm));

	% compute LFP and firing rate lags
	for i = 1:length(paradigm)
		S = PID(a:z,i);
		X = -raw_LFP(a:z,i);
		R = fA(a:z,i);
		X_pred = LFP_model(a:z,i);

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

		% compute stimulus auto-correlation 
		temp = mean(S_acorr,2);
		temp = temp(1e3:end);
		try
			S_autocorr(i) = find(temp<1/exp(1),1,'first');
		catch
		end

	end

	if di == 1
		ax(1) = subplot(2,3,1); hold on
		ax(2) = subplot(2,3,2); hold on
		% compare LFP cross correlations at the lowest and highest dose
		c = parula(max(paradigm)+1);
		temp = LFP_xcorr(:,:,paradigm==1);
		temp = reshape(temp,2*chunk_size - 1, size(temp,2)*size(temp,3));
		lags = 1e-3*(1:length(temp)) - 1;
		plot(ax(1),lags,mean(temp,2)/max(mean(temp,2)),'Color',c(1,:))
		temp = LFP_xcorr(:,:,paradigm==10);
		temp = reshape(temp,2*chunk_size - 1, size(temp,2)*size(temp,3));
		lags = 1e-3*(1:length(temp)) - 1;
		plot(ax(1),lags,mean(temp,2)/max(mean(temp,2)),'Color',c(10,:))


		% compare firing cross correlations at the lowest and highest dose
		c = parula(max(paradigm)+1);
		temp = fA_xcorr(:,:,paradigm==1);
		temp = reshape(temp,2*chunk_size - 1, size(temp,2)*size(temp,3));
		lags = 1e-3*(1:length(temp)) - 1;
		plot(ax(2),lags,mean(temp,2)/max(mean(temp,2)),'Color',c(1,:))
		temp = fA_xcorr(:,:,paradigm==10);
		temp = reshape(temp,2*chunk_size - 1, size(temp,2)*size(temp,3));
		lags = 1e-3*(1:length(temp)) - 1;
		plot(ax(2),lags,nanmean(temp,2)/max(nanmean(temp,2)),'Color',c(10,:))



	end


	firing_lags(firing_lags<-100) = NaN; % obviously wrong


	ax(2+di) = subplot(2,3,2+di); hold on

	% show lags of firing rate and LFP for all doses
	mean_stim = mean(PID(a:z,:));
	if di == 1
		LFP_lags(mean_stim > .7 & mean_stim < .9) = NaN; % firingr rates for this data don't exist
	end
	c = lines(10);
	LFP_color = c(5,:);
	firing_color = c(4,:);
	for i = 1:max(paradigm)
		x = nonnans(firing_lags(paradigm==i));
		errorbar(ax(di+2),mean(mean_stim(paradigm==i)),mean(x),sem(x),'Color',LFP_color)
		x = nonnans(LFP_lags(paradigm==i));
		errorbar(ax(di+2),mean(mean_stim(paradigm==i)),mean(x),sem(x),'Color',firing_color)

		x = nonnans(LFP_model_lags(paradigm==i));
		errorbar(ax(di+2),mean(mean_stim(paradigm==i)),mean(x),sem(x),'Color','r')

	end
	%plot(mean_stim,firing_lags,'.','Color',firing_color,'MarkerSize',20);
	%plot(mean_stim,LFP_lags,'.','Color',LFP_color,'MarkerSize',20);
	xlabel(ax(2+di),'\mu_{Stimulus} (V)')
	ylabel(ax(2+di),'Lag (ms)')


	[LFP_model_lags,idx] = sort(LFP_model_lags);
	plot(mean_stim(idx),LFP_model_lags,'r.-')


	legend({'Firing rate','LFP','LFP model'},'Location','northwest')
	set(ax(2+di),'XLim',[0 x_limits(di)],'YLim',[0 200])
	title([orn_names{di} ' - ' odour_names{di}])

end

set(ax(1),'XLim',[0 1],'YLim',[-.5 1.1])
set(ax(2),'XLim',[0 1],'YLim',[-.5 1.1])
xlabel(ax(1),'Lag (s)')
xlabel(ax(2),'Lag (s)')
ylabel(ax(1),'Cross correlation (norm)')
ylabel(ax(2),'Cross correlation (norm)')
title(ax(1),'Stimulus \rightarrow LFP')
title(ax(2),'Stimulus \rightarrow Firing rate')


prettyFig;

if being_published	
	snapnow	
	delete(gcf)
end


%% Version Info
%
pFooter;




