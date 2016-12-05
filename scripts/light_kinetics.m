% fig_light_mechanism.m
% 
% created by Srinivas Gorur-Shandilya at 9:50 , 28 December 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


%% Kinetics of response with light or odor background
% In the following figure, I plot the kinetics of response as a function of light or odor background. All data comes from Chrimson expressing ab3A ORNs. On the left, I plot the LFP and firing rate lags w.r.t a fluctuating odor stimulus as I increase the background light intensity. In the panel on the right, I plot the firing rate lag w.r.t to a fluctuating light stimulus as a function of the background odor concentration. No clear trend is visible anywhere. 

pHeader;
dm = dataManager;

% build the LED map
x = [.75:0.05:1.1 1.2:.1:3.6 3.8 4 4.2 4.5];
x = [0:0.1:0.7 x];
y = [0 4 12.8 24.1 36.2 48.1 60 70 95 114 134 151 167 184 201 219 236 252 269 283 299 314 328 343 357 370 383 394 404 413 419 422 422 421 419 418 417];
y = [0*(0:0.1:0.7) y ];
light_power_fit = fit(x(:),y(:),'spline');



			;;       ;;;;  ;;;;;;   ;;     ;; ;;;;;;;;    ;;;;;;;;   ;;;;;;   
			;;        ;;  ;;    ;;  ;;     ;;    ;;       ;;     ;; ;;    ;;  
			;;        ;;  ;;        ;;     ;;    ;;       ;;     ;; ;;        
			;;        ;;  ;;   ;;;; ;;;;;;;;;    ;;       ;;;;;;;;  ;;   ;;;; 
			;;        ;;  ;;    ;;  ;;     ;;    ;;       ;;     ;; ;;    ;;  
			;;        ;;  ;;    ;;  ;;     ;;    ;;       ;;     ;; ;;    ;;  
			;;;;;;;; ;;;;  ;;;;;;   ;;     ;;    ;;       ;;;;;;;;   ;;;;;;   
			 

clearvars -except ax axs being_published dm light_power_fit od od_odour uts main_fig supp_fig

[PID, LFP, fA, paradigm, orn, fly, AllControlParadigms, paradigm_hashes] = consolidateData(dm.getPath('9df40545ca5197f63f2dd4af7c905316'),true);
a = 15e3; z = 55e3;

% clean up the data
fA(:,sum(fA) == 0) = NaN;

% figure out which paradigms to plot
light_background_paradigms = false(length(paradigm),1);
odour_background_paradigms = false(length(paradigm),1);
for i = 1:length(paradigm)
	if any(strfind(AllControlParadigms(paradigm(i)).Name,'0%')) || any(strfind(AllControlParadigms(paradigm(i)).Name,'Light'))
		light_background_paradigms(i) = true;
	end
	if any(strfind(AllControlParadigms(paradigm(i)).Name,'%')) 
		odour_background_paradigms(i) = true;
	end
end

% calculate the mean light power in every trial
mean_light_power = zeros(length(paradigm),1);
mean_PID = mean(PID(a:z,:));
for i = 1:length(paradigm)
	this_p = paradigm(i);
	try
		mean_light_power(i) = light_power_fit(str2double(AllControlParadigms(this_p).Name(7:strfind(AllControlParadigms(this_p).Name,'V')-1)));
	catch
	end
end
mean_light_power(isnan(mean_light_power)) = 0;


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
	if light_background_paradigms(i)
		S = PID(a:z,i);
		X = -LFP(a:z,i);
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
end


firing_lags(firing_lags<0) = NaN; % obviously wrong
LFP_lags(LFP_lags<0) = NaN; % obviously wrong

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1); hold on
% show lags of firing rate and LFP for all doses
mean_stim = mean_light_power;
for i = 1:max(paradigm)
	x = nonnans(firing_lags(paradigm==i));
	errorbar(mean(mean_stim(paradigm==i)),mean(x),sem(x),'k')
	x = nonnans(LFP_lags(paradigm==i));
	errorbar(mean(mean_stim(paradigm==i)),mean(x),sem(x),'r')
end
xlabel('\mu_{Stimulus} (\muW)')
ylabel('Lag (ms)')
legend({'Firing rate','LFP'},'Location','northwest')
set(gca,'XLim',[-50 450],'YLim',[0 200])



			 ;;;;;;;  ;;;;;;;;   ;;;;;;;  ;;;;;;;;     ;;;;;;;;   ;;;;;;   
			;;     ;; ;;     ;; ;;     ;; ;;     ;;    ;;     ;; ;;    ;;  
			;;     ;; ;;     ;; ;;     ;; ;;     ;;    ;;     ;; ;;        
			;;     ;; ;;     ;; ;;     ;; ;;;;;;;;     ;;;;;;;;  ;;   ;;;; 
			;;     ;; ;;     ;; ;;     ;; ;;   ;;      ;;     ;; ;;    ;;  
			;;     ;; ;;     ;; ;;     ;; ;;    ;;     ;;     ;; ;;    ;;  
			 ;;;;;;;  ;;;;;;;;   ;;;;;;;  ;;     ;;    ;;;;;;;;   ;;;;;;   

if ~exist('od','var')
	p = getPath(dataManager,'6c02fb233f8221b2d22d98ebf8a48edb');
	od = raw2ORNData(p,'led_power_func',light_power_fit,'use_led',true,'filter_LFP',false);

	% also create another object with the odour values
	od_odour = raw2ORNData(p,'led_power_func',light_power_fit,'use_led',false,'filter_LFP',false);
end

% reshape
LED = [od.stimulus];
PID = [od_odour.stimulus] - min([od_odour(:,1).stimulus]);
fA = [od.firing_rate];
paradigm = 0*(1:size(LED,2));

firing_lags = NaN*paradigm;
firing_max_corr = NaN*paradigm;
a = 35e3+1; z = 55e3;
chunk_size = 1e3;

fA_xcorr = NaN(chunk_size*2 - 1,20,length(paradigm));

% compute LFP and firing rate lags
for i = 1:length(paradigm)
	S = LED(a:z,i);
	R = fA(a:z,i);

	% reshape into chunks
	S = reshape(S,chunk_size,length(S)/chunk_size);
	R = reshape(R,chunk_size,length(R)/chunk_size);

	S = bsxfun(@minus, S, mean(S));
	R = bsxfun(@minus, R, mean(R));

	R_lag = NaN(chunk_size*2-1,size(S,2));
	clear R_lag
	for j = 1:size(S,2)
		R_lag(:,j) = xcorr(R(:,j)/std(R(:,j)),S(:,j)/std(S(:,j)));
	end

	R_lag = R_lag/chunk_size;
	fA_xcorr(:,:,i) = R_lag;
	R_lag = mean(R_lag,2);
	[firing_max_corr(i),loc] = max(R_lag);
	firing_lags(i) = loc - 1e3;
end

subplot(1,2,2); hold on
% show lags of firing rate and LFP for all doses
mean_stim = mean(PID(a:z,:));
plot(mean_stim,firing_lags,'k+')
xlabel('\mu_{Stimulus} (V)')
ylabel('Firing Lag (ms)')
set(gca,'XLim',[0 .30],'YLim',[0 200])

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%% Version Info
%
pFooter;


%% Supplementary Figure
