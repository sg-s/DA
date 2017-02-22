


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


data_hashes = {'93ba5d68174e3df9f462a1fc48c581da','cd6753c0e4cf02895cd5e2c5cb58aa1a','3ea08ccfa892c6545d74bbdaaa6cbee1'};
odour_names = {'ethyl-acetate','1-pentanol','2-butanone'};
orn_names = {'ab3A','ab2A','ab2A'};
x_limits = [1.7 .15 .85 ];

for di = 1:length(data_hashes)

	MSGdata = consolidateData2(getPath(dataManager,data_hashes{di}));
	MSGdata = cleanMSGdata(MSGdata,'extract_filter',false);
	v2struct(MSGdata)

	% remove baseline from PID for each trial
	PID = bsxfun(@minus,PID,min(PID));


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

		% reshape into chunks
		S = reshape(S,chunk_size,length(S)/chunk_size);
		X = reshape(X,chunk_size,length(X)/chunk_size);
		R = reshape(R,chunk_size,length(R)/chunk_size);

		S = bsxfun(@minus, S, mean(S));
		X = bsxfun(@minus, X, mean(X));
		R = bsxfun(@minus, R, mean(R));

		X_lag = NaN(chunk_size*2-1,size(S,2));
		R_lag = NaN(chunk_size*2-1,size(S,2));
		S_acorr = NaN(chunk_size*2-1,size(S,2));
		for j = 1:size(S,2)
			X_lag(:,j) = xcorr(X(:,j)/std(X(:,j)),S(:,j)/std(S(:,j)));
			R_lag(:,j) = xcorr(R(:,j)/std(R(:,j)),S(:,j)/std(S(:,j)));
			S_acorr(:,j) = xcorr(S(:,j)/std(S(:,j)),S(:,j)/std(S(:,j)));
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
	
	c = lines(10);
	LFP_color = c(5,:);
	firing_color = c(4,:);
	
	if di == 1
		x = NaN*(1:max(paradigm));
		l = NaN*(1:max(paradigm));
		l_e = NaN*(1:max(paradigm));
		f = NaN*(1:max(paradigm));
		f_e = NaN*(1:max(paradigm));
		for i = 1:max(paradigm)
			temp = nonnans(firing_lags(paradigm==i));
			f(i) = mean(temp);
			f_e(i) = sem(temp);
			x(i) = mean(mean_stim(paradigm==i));
			temp = nonnans(LFP_lags(paradigm==i));
			l(i) = mean(temp);
			l_e(i) = sem(temp);
		end
		errorbar(x,l,l_e,'Color',LFP_color);
		errorbar(x,f,f_e,'Color',firing_color);
	else
		plotPieceWiseLinear(mean_stim,LFP_lags,'Color',LFP_color,'nbins',10);
		plotPieceWiseLinear(mean_stim,firing_lags,'Color',firing_color,'nbins',10);
	end

	
	xlabel(ax(2+di),'\mu_{Stimulus} (V)')
	ylabel(ax(2+di),'Lag (ms)')

	legend({'LFP','Firing rate'},'Location','northwest')
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


% now show the same thing in the naturalistic stimulus date 
min_acceptable_corr = .5;
min_acceptable_lag = 10;

subplot(2,3,6); hold on
load(getPath(dataManager,'aeb361c027b71938021c12a6a12a85cd'),'-mat')

PID = nanmean([od.stimulus],2);
LFP = nanmean([od.LFP],2);
fA = [od.firing_rate];
fA(:,sum(fA)==0) = [];
fA = nanmean(fA,2);

[lag, mean_x, max_corr] = findLagAndMeanInWindow(PID(5e3:end-5e3),-LFP(5e3:end-5e3),1e3,50);
rm_this = lag<min_acceptable_lag | max_corr < min_acceptable_corr;
lag(rm_this) = [];
mean_x(rm_this) = [];
plotPieceWiseLinear(mean_x,lag,'Color',LFP_color,'nbins',10);

[lag, mean_x, max_corr] = findLagAndMeanInWindow(PID(5e3:end-5e3),fA(5e3:end-5e3),1e3,50);
rm_this = lag<min_acceptable_lag | max_corr < min_acceptable_corr;
lag(rm_this) = [];
mean_x(rm_this) = [];
plotPieceWiseLinear(mean_x,lag,'Color',firing_color,'nbins',10);

xlabel('\mu_{Stimulus} in preceding 1s')
ylabel('Lag (ms)')
set(gca,'YLim',[0 200])
title(['ab3A - ethyl acetate' char(10) 'naturalistic stimulus'])

prettyFig;

labelFigure('x_offset',.01)


if being_published	
	snapnow	
	delete(gcf)
end


%% Version Info
%
pFooter;




