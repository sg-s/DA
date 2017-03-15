


pHeader;

%% Kinetics of LFP and firing rate 
% In this document, I analyse the kinetics of firing rate and LFP w.r.t the stimulus using cross correlation functions. 

%%
% In the following figure, I how the lags of the LFP and firing rate vary with stimulus background from the gaussian fluctuation data. (a) Normalised cross correlation computed from stimulus to LFP for the lowest dose (purple) and for the highest dose (yellow). (b) Raw LFP traces at odour offset show that LFP kinetics slow down with increasing mean stimulus (brighter colours). (c) Normalised cross correlation computed from stimulus to firing rate for the lowest dose (purple) and for the highest dose (yellow). No major change is visible. (d) Lags of firing rate (black) and LFP (red) estimated from peaks of cross correlation functions, as a function of the mean stimulus in that trial. Error bars are standard error of the mean. Cross correlations were estimated from 20 seconds of data in each trial, chunked into 1-second long fragments, and then averaged. 

c = lines(5);
LFP_color = c(4,:);
firing_color = c(5,:);

main_fig = figure('outerposition',[0 0 1400 801],'PaperUnits','points','PaperSize',[1400 801]); hold on
ax.lfp_xcorr = subplot(2,4,1); hold on
ax.firing_xcorr = subplot(2,4,2); hold on
ax.autocorr1 = subplot(2,4,3); hold on
ax.odor_lags(1) = subplot(2,4,4); hold on
ax.odor_lags(2) = subplot(2,4,5); hold on
ax.odor_lags(3) = subplot(2,4,6); hold on
ax.odor_lags(4) = subplot(2,4,7); hold on
ax.light_lags = subplot(2,4,8); hold on

% data_hashes = {'93ba5d68174e3df9f462a1fc48c581da','cd6753c0e4cf02895cd5e2c5cb58aa1a','3ea08ccfa892c6545d74bbdaaa6cbee1'};
% odour_names = {'ethyl-acetate','1-pentanol','2-butanone'};
% orn_names = {'ab3A','ab2A','ab2A'};
% x_limits = [1.7 .15 .85 ];


% define what we want to work on

% define what we want to work on
data_hashes = {'93ba5d68174e3df9f462a1fc48c581da','bcd4cf4fe12817d084a2b06f981161ee','cd6753c0e4cf02895cd5e2c5cb58aa1a','3ea08ccfa892c6545d74bbdaaa6cbee1'};
odour_names = {'ethyl-acetate','1-pentanol','1-pentanol','2-butanone'};
orn_names = {'ab3A','ab3A','ab2A','ab2A'};
x_limits = [1.7 .035 .2 .8];
nbins = [10 5 10 10 10];

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

		% compare LFP cross correlations at the lowest and highest dose
		c = parula(max(paradigm)+1);
		temp = LFP_xcorr(:,:,paradigm==1);
		temp = reshape(temp,2*chunk_size - 1, size(temp,2)*size(temp,3));
		lags = 1e-3*(1:length(temp)) - 1;
		plot(ax.lfp_xcorr,lags,mean(temp,2)/max(mean(temp,2)),'Color',c(1,:))
		temp = LFP_xcorr(:,:,paradigm==10);
		temp = reshape(temp,2*chunk_size - 1, size(temp,2)*size(temp,3));
		lags = 1e-3*(1:length(temp)) - 1;
		plot(ax.lfp_xcorr,lags,mean(temp,2)/max(mean(temp,2)),'Color',c(10,:))


		% compare firing cross correlations at the lowest and highest dose
		c = parula(max(paradigm)+1);
		temp = fA_xcorr(:,:,paradigm==1);
		temp = reshape(temp,2*chunk_size - 1, size(temp,2)*size(temp,3));
		lags = 1e-3*(1:length(temp)) - 1;
		plot(ax.firing_xcorr,lags,mean(temp,2)/max(mean(temp,2)),'Color',c(1,:))
		temp = fA_xcorr(:,:,paradigm==10);
		temp = reshape(temp,2*chunk_size - 1, size(temp,2)*size(temp,3));
		lags = 1e-3*(1:length(temp)) - 1;
		plot(ax.firing_xcorr,lags,nanmean(temp,2)/max(nanmean(temp,2)),'Color',c(10,:))

		% also compare auto-correlations of the stimulus at high and low 
		c = parula(max(paradigm)+1);
		temp = S_a(:,:,paradigm==1);
		temp = reshape(temp,2*chunk_size - 1, size(temp,2)*size(temp,3));
		lags = 1e-3*(1:length(temp)) - 1;
		plot(ax.autocorr1,lags,mean(temp,2)/max(mean(temp,2)),'Color',c(1,:))
		temp = S_a(:,:,paradigm==10);
		temp = reshape(temp,2*chunk_size - 1, size(temp,2)*size(temp,3));
		lags = 1e-3*(1:length(temp)) - 1;
		plot(ax.autocorr1,lags,nanmean(temp,2)/max(nanmean(temp,2)),'Color',c(10,:))


	end


	firing_lags(firing_lags<-100) = NaN; % obviously wrong


	% show lags of firing rate and LFP for all doses
	mean_stim = mean(PID(a:z,:));
	
	c = lines(10);
	LFP_color = c(5,:);
	firing_color = c(4,:);
	
	[lfp_stats(di).rho,lfp_stats(di).p] = spear(mean_stim,LFP_lags);
	[firing_stats(di).rho,firing_stats(di).p] = spear(mean_stim,firing_lags);

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
		errorbar(ax.odor_lags(di),x,l,l_e,'Color',LFP_color);
		errorbar(ax.odor_lags(di),x,f,f_e,'Color',firing_color);
	else
		% bad_trials = LFP_lags < 10;
		% LFP_lags(bad_trials) = NaN;
		% firing_lags(bad_trials) = NaN;
		plotPieceWiseLinear(ax.odor_lags(di),mean_stim,LFP_lags,'Color',LFP_color,'nbins',nbins(di));
		plotPieceWiseLinear(ax.odor_lags(di),mean_stim,firing_lags,'Color',firing_color,'nbins',nbins(di));
	end

	xlabel(ax.odor_lags(di),'\mu_{Stimulus} (V)')
	ylabel(ax.odor_lags(di),'Lag (ms)')

	set(ax.odor_lags(di),'XLim',[0 x_limits(di)],'YLim',[0 200])
	title(ax.odor_lags(di),[orn_names{di} ' - ' odour_names{di}])

end

clear l
l(1) = plot(ax.odor_lags(1),NaN,NaN,'Color',LFP_color,'Marker','.','MarkerSize',24,'LineStyle','none');
l(2) = plot(ax.odor_lags(1),NaN,NaN,'Color',firing_color,'Marker','.','MarkerSize',24,'LineStyle','none');
legend(l,{'LFP','Firing rate'},'Location','northwest')


set(ax.lfp_xcorr,'XLim',[0 1],'YLim',[-.5 1.1])
set(ax.firing_xcorr,'XLim',[0 1],'YLim',[-.5 1.1])
set(ax.autocorr1,'XLim',[0 1],'YLim',[-.5 1.1])
xlabel(ax.lfp_xcorr,'Lag (s)')
xlabel(ax.firing_xcorr,'Lag (s)')
xlabel(ax.autocorr1,'Lag (s)')
ylabel(ax.lfp_xcorr,'Cross correlation (norm)')
ylabel(ax.firing_xcorr,'Cross correlation (norm)')
ylabel(ax.autocorr1,'Autocorrelation (norm)')
title(ax.lfp_xcorr,'Stimulus \rightarrow LFP')
title(ax.firing_xcorr,'Stimulus \rightarrow Firing rate')
title(ax.autocorr1,'Stimulus auto-correlation')

% add some statistics for each plot 
for i = 1:4
	rho = lfp_stats(i).rho;
	p = lfp_stats(i).p;
	if p < 1e-4 
		tp = 'p<10^-^4';
	else
		tp = ['p = ' oval(p)];
	end
	t = ['\rho = ' oval(rho), ', ' tp];
	th(i) = text(.4,100,t);
	th(i).Parent = ax.odor_lags(i);
	th(i).Color = LFP_color;
	th(i).FontSize = 18;
end

th(1).Position = [.2 142];
th(2).Position = [.005 160];
th(3).Position = [.05 160];
th(4).Position = [.1 140];
% th(5).Position = [.05 160];

% add some statistics for each plot 
for i = 1:4
	rho = firing_stats(i).rho;
	p = firing_stats(i).p;
	if p < 1e-4 
		tp = 'p<10^-^4';
	else
		tp = ['p = ' oval(p)];
	end
	t = ['\rho = ' oval(rho), ', ' tp];
	th2(i) = text(.4,100,t);
	th2(i).Parent = ax.odor_lags(i);
	th2(i).Color = firing_color;
	th2(i).FontSize = 18;
end

th2(1).Position = [.2 20];
th2(2).Position = [.005 40];
th2(3).Position = [.05 20];
th2(4).Position = [.1 20];
% th2(5).Position = [.05 20];


% now show that it speeds up with light 

;;       ;;;;  ;;;;;;   ;;     ;; ;;;;;;;; 
;;        ;;  ;;    ;;  ;;     ;;    ;;    
;;        ;;  ;;        ;;     ;;    ;;    
;;        ;;  ;;   ;;;; ;;;;;;;;;    ;;    
;;        ;;  ;;    ;;  ;;     ;;    ;;    
;;        ;;  ;;    ;;  ;;     ;;    ;;    
;;;;;;;; ;;;;  ;;;;;;   ;;     ;;    ;;    


dm = dataManager;

% build the LED map
x = [.75:0.05:1.1 1.2:.1:3.6 3.8 4 4.2 4.5];
x = [0:0.1:0.7 x];
y = [0 4 12.8 24.1 36.2 48.1 60 70 95 114 134 151 167 184 201 219 236 252 269 283 299 314 328 343 357 370 383 394 404 413 419 422 422 421 419 418 417];
y = [0*(0:0.1:0.7) y ];
light_power_fit = fit(x(:),y(:),'spline');

[~, ~, fA, paradigm, orn, fly, AllControlParadigms] = consolidateData(dm.getPath('76ad69b63717ab07183afe0bba30cb4d'),true);


% make new vectors for mean and range of the stimulus
a = 25e3; z = 55e3;
LED = NaN*fA;
r = NaN*paradigm; m = NaN*paradigm;
light_s = NaN*paradigm; light_m = NaN*paradigm;
for i = 1:length(paradigm)
	temp = AllControlParadigms(paradigm(i)).Name;
	r(i) = (str2double(temp(3:strfind(temp,'_')-1)));
	m(i) = (str2double(temp(strfind(temp,'_')+3:end)));
	LED(:,i) = light_power_fit(AllControlParadigms(paradigm(i)).Outputs(1,1:10:end));
	light_s(i) = std(LED(a:z,i));
	light_m(i) = mean(LED(a:z,i));
end

% find delays for every trial
lags = NaN*m;
warning off
for i = 1:length(lags)
	lags(i) = finddelay(LED(a:z,i),fA(a:z,i));
end
warning on
lags(lags==0) = NaN;


ll = unique(light_m);
x = NaN*ll; y = NaN*ll; ye = NaN*ll;
for i = 1:length(ll)
	this = light_m == ll(i);
	x(i) = ll(i);
	y(i) = nanmean(lags(this));
	ye(i) = nanstd(lags(this))/sqrt(sum(this));
end
errorbar(ax.light_lags,x,y,ye,'Color','k')

xlabel(ax.light_lags,'Mean light power (\muW)')
ylabel(ax.light_lags,'Firing lag (ms)')
set(ax.light_lags,'YLim',[55 75],'XLim',[0 150])
title(ax.light_lags,'ab3A - light')

% add some statistics for each plot 
[rho,p] = spear(light_m,lags);
tp = 'p < 10^-^4';
t = ['\rho = ' oval(rho), ', ' tp];
th3 = text(.4,100,t);
th3.Parent = ax.light_lags;
th3.Color = [0 0 0];
th3.FontSize = 18;
th3.Position = [50 58];

% add legend to the first plot showing what the purple and yellow curves are
c = parula(11);
axes(ax.lfp_xcorr)
th1 = text(.7,1,'Low stimulus','Color',c(1,:));
th2 = text(.7,.8,'High stimulus','Color',c(10,:));

th1.Position = [.5 1];
th2.Position = [.5 .88];

figure(main_fig)
prettyFig;

labelFigure('x_offset',-.01)


if being_published	
	snapnow	
	delete(gcf)
end




% dm = dataManager;

% % build the LED map
% x = [.75:0.05:1.1 1.2:.1:3.6 3.8 4 4.2 4.5];
% x = [0:0.1:0.7 x];
% y = [0 4 12.8 24.1 36.2 48.1 60 70 95 114 134 151 167 184 201 219 236 252 269 283 299 314 328 343 357 370 383 394 404 413 419 422 422 421 419 418 417];
% y = [0*(0:0.1:0.7) y ];
% light_power_fit = fit(x(:),y(:),'spline');

% [~, ~, fA, paradigm, orn, fly, AllControlParadigms] = consolidateData(dm.getPath('38901017007c50ea52637891619ab91c'),true);


% % make new vectors for mean and range of the stimulus
% a = 25e3; z = 55e3;
% LED = NaN*fA;
% r = NaN*paradigm; m = NaN*paradigm;
% light_s = NaN*paradigm; light_m = NaN*paradigm;
% for i = 1:length(paradigm)
% 	temp = AllControlParadigms(paradigm(i)).Name;
% 	r(i) = (str2double(temp(3:strfind(temp,'_')-1)));
% 	m(i) = (str2double(temp(strfind(temp,'_')+3:end)));
% 	LED(:,i) = light_power_fit(AllControlParadigms(paradigm(i)).Outputs(1,1:10:end));
% 	light_s(i) = std(LED(a:z,i));
% 	light_m(i) = mean(LED(a:z,i));
% end

% % find delays for every trial
% lags = NaN*m;
% warning off
% for i = 1:length(lags)
% 	lags(i) = finddelay(LED(a:z,i),fA(a:z,i));
% end
% warning on
% lags(lags==0) = NaN;

% plotPieceWiseLinear(ax.light_lags,light_m,lags,'nbins',5);
% xlabel(ax.light_lags,'Mean light power (\muW)')
% ylabel(ax.light_lags,'Firing lag (ms)')
% set(ax.light_lags,'YLim',[60 80])
% title(ax.light_lags,'ab3A - light')

% % add some statistics for each plot 
% [rho,p] = spear(light_m,lags);
% tp = 'p < 10^-^4';
% t = ['\rho = ' oval(rho), ', ' tp];
% th3 = text(.4,100,t);
% th3.Parent = ax.light_lags;
% th3.Color = [0 0 0];
% th3.FontSize = 18;
% th3.Position = [100 75];



%% Version Info
%
pFooter;




