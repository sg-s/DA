% Calcium.m
% 
% created by Srinivas Gorur-Shandilya at 4:00 , 04 November 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

%% Calcium Manipulations of Extracellular Medium
% In this document we manipulate the extra-cellular calcium concentration using EDTA (1mM, to lower it), or increasing the Calcium concentration to 10mM (from a normal conc. of 1mM). 



pHeader;

[alldata(1).PID, alldata(1).LFP, alldata(1).fA, alldata(1).paradigm, alldata(1).orn, alldata(1).fly, alldata(1).AllControlParadigms, alldata(1).paradigm_hashes] = consolidateData('/Volumes/sgs_data/Dropbox (emonetlab)/users/srinivas_gs/data/DA-project/electrode-calcium/low/',true);

[alldata(2).PID, alldata(2).LFP, alldata(2).fA, alldata(2).paradigm, alldata(2).orn, alldata(2).fly, alldata(2).AllControlParadigms, alldata(2).paradigm_hashes] = consolidateData('/Volumes/sgs_data/Dropbox (emonetlab)/users/srinivas_gs/data/DA-project/electrode-calcium/normal/',true);

[alldata(3).PID, alldata(3).LFP, alldata(3).fA, alldata(3).paradigm, alldata(3).orn, alldata(3).fly, alldata(3).AllControlParadigms, alldata(3).paradigm_hashes] = consolidateData('/Volumes/sgs_data/Dropbox (emonetlab)/users/srinivas_gs/data/DA-project/electrode-calcium/high/',true);

for ai = 1:length(alldata)
	% remove baseline from all PIDs
	for i = 1:width(alldata(ai).PID)
		alldata(ai).PID(:,i) = alldata(ai).PID(:,i) - mean(alldata(ai).PID(1:5e3,i));
	end

	% remove baseline from all LFPs
	for i = 1:width(alldata(ai).LFP)
		alldata(ai).LFP(:,i) = alldata(ai).LFP(:,i) - mean(alldata(ai).LFP(1:5e3,i));
	end

	% band pass all the LFP
	alldata(ai).filtered_LFP = alldata(ai).LFP;
	for i = 1:width(alldata(ai).LFP)
		alldata(ai).filtered_LFP(:,i) = bandPass(alldata(ai).LFP(:,i),1000,10);
	end

	% remove "Flicker" from paradigm names
	for i = 1:length(alldata(ai).AllControlParadigms)
		alldata(ai).AllControlParadigms(i).Name = strrep(alldata(ai).AllControlParadigms(i).Name,'Flicker-','');
	end

end

%% Baseline Firing
% In this section we look at the baseline firing of the neuron. 

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
for ai = 1:length(alldata)
	temp = nonnans(nonzeros(mean(alldata(ai).fA(1:5e3,:))));
	errorbar(ai,mean(temp),sem(temp),'k')
end
set(gca,'XTick',1:length(alldata),'XTickLabel',{'Low Ca^{2+}','Normal Ca^{2+}','High Ca^{2+}'},'XTickLabelRotation',45)
ylabel('Baseline Firing Rate (Hz)')
prettyFig()

if being_published
	snapnow
	delete(gcf)
end


%% LFP Gain
% In this section we compare the LFP response to a fluctuating stimulus between normal calcium and low and high manipulations. 

a = 10e3;
z = 45e3;
ss = 20;
figure('outerposition',[0 0 1200 700],'PaperUnits','points','PaperSize',[1200 700]); hold on
for i = 1:3
	subplot(2,3,i), hold on
	normal_LFP =  alldata(2).filtered_LFP(a:z,alldata(2).paradigm==i);
	test_LFP = alldata(1).filtered_LFP(a:z,alldata(1).paradigm==i);
	m = 1.2*min([test_LFP(:); normal_LFP(:)]);
	M = 1.2*max([test_LFP(:); normal_LFP(:)]);
	plot([m M],[m M],'k--')
	plot(normal_LFP(1:ss:end),test_LFP(1:ss:end),'b.')
	axis square
	xlabel('\Delta LFP at normal Ca^{2+}')
	ylabel('\Delta LFP at low Ca^{2+}')
	title(alldata(1).AllControlParadigms(find(alldata(2).paradigm==i)).Name)

	subplot(2,3,3+i), hold on
	test_LFP = alldata(3).filtered_LFP(a:z,alldata(3).paradigm==i);
	m = 1.2*min([test_LFP(:); normal_LFP(:)]);
	M = 1.2*max([test_LFP(:); normal_LFP(:)]);
	plot([m M],[m M],'k--')
	plot(normal_LFP(1:ss:end),test_LFP(1:ss:end),'r.')
	axis square
	xlabel('\Delta LFP at normal Ca^{2+}')
	ylabel('\Delta LFP at high Ca^{2+}')
end

prettyFig()

if being_published
	snapnow
	delete(gcf)
end

%% Firing Gain 
% In this section we compare the Firing response to a fluctuating stimulus between normal calcium and low and high manipulations. 

a = 10e3;
z = 45e3;
ss = 20;
figure('outerposition',[0 0 1200 700],'PaperUnits','points','PaperSize',[1200 700]); hold on
for i = 1:3
	subplot(2,3,i), hold on
	normal_LFP =  alldata(2).fA(a:z,alldata(2).paradigm==i);
	test_LFP = alldata(1).fA(a:z,alldata(1).paradigm==i);
	m = 0;
	M = 100;
	plot([m M],[m M],'k--')
	plot(normal_LFP(1:ss:end),test_LFP(1:ss:end),'b.')
	axis square
	xlabel('Firing Rate at normal Ca^{2+}')
	ylabel('Firing Rate at low Ca^{2+}')
	title(alldata(1).AllControlParadigms(find(alldata(2).paradigm==i)).Name)

	subplot(2,3,3+i), hold on
	test_LFP = alldata(3).fA(a:z,alldata(3).paradigm==i);
	m = 0;
	M = 100;
	plot([m M],[m M],'k--')
	plot(normal_LFP(1:ss:end),test_LFP(1:ss:end),'r.')
	axis square
	xlabel('Firing Rate at normal Ca^{2+}')
	ylabel('Firing Rate at high Ca^{2+}')
end

prettyFig()

if being_published
	snapnow
	delete(gcf)
end

;;       ;;;;;;;; ;;;;;;;;     ;;    ;; ;;;; ;;    ;; ;;;;;;;; ;;;;;;;; ;;;;  ;;;;;;   ;;;;;;  
;;       ;;       ;;     ;;    ;;   ;;   ;;  ;;;   ;; ;;          ;;     ;;  ;;    ;; ;;    ;; 
;;       ;;       ;;     ;;    ;;  ;;    ;;  ;;;;  ;; ;;          ;;     ;;  ;;       ;;       
;;       ;;;;;;   ;;;;;;;;     ;;;;;     ;;  ;; ;; ;; ;;;;;;      ;;     ;;  ;;        ;;;;;;  
;;       ;;       ;;           ;;  ;;    ;;  ;;  ;;;; ;;          ;;     ;;  ;;             ;; 
;;       ;;       ;;           ;;   ;;   ;;  ;;   ;;; ;;          ;;     ;;  ;;    ;; ;;    ;; 
;;;;;;;; ;;       ;;           ;;    ;; ;;;; ;;    ;; ;;;;;;;;    ;;    ;;;;  ;;;;;;   ;;;;;;  


%% Kinetics
% In this section, I look at how the kinetics of the LFP changes as a function of the extracellular calcium. 

a = 35e3+1; z = 55e3;
chunk_size = 1e3;

for i = 1:length(alldata)
	paradigm = alldata(i).paradigm;
	alldata(i).LFP_lags = NaN*paradigm;
	alldata(i).LFP_max_corr = NaN*paradigm;
	alldata(i).firing_lags = NaN*paradigm;
	alldata(i).firing_max_corr = NaN*paradigm;
	alldata(i).LFP_xcorr = NaN(chunk_size*2 - 1,20,length(paradigm));
	alldata(i).fA_xcorr = NaN(chunk_size*2 - 1,20,length(paradigm));
end

% compute LFP and firing rate lags
for k = 1:length(alldata)

	PID = alldata(k).PID;
	LFP = alldata(k).LFP;
	fA = alldata(k).fA;
	paradigm = alldata(k).paradigm;

	for i = 1:length(paradigm)

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
		alldata(k).LFP_xcorr(:,:,i) = X_lag;
		X_lag = mean(X_lag,2);
		[alldata(k).LFP_max_corr(i),loc] = max(X_lag);
		alldata(k).LFP_lags(i) = loc - 1e3;

		R_lag = R_lag/chunk_size;
		alldata(k).fA_xcorr(:,:,i) = R_lag;
		R_lag = mean(R_lag,2);
		[alldata(k).firing_max_corr(i),loc] = max(R_lag);
		alldata(k).firing_lags(i) = loc - 1e3;
	end
end

figure('outerposition',[0 0 1501 802],'PaperUnits','points','PaperSize',[1501 802]); hold on
c = [0 0 1; 0 0 0; 1 0 0];
for i = 1:3 % because we have 3 values
	subplot(2,3,i); hold on
	set(gca,'XLim',[-100 500])
	xlabel('Lag (ms)')
	ylabel('Cross correlation (norm)')
	title(['\mu_{Stimulus} = ' oval(mean(alldata(1).PID(a:z,i)))])
	% for each calcium conc, plot the xcorr b/w stimulus and LFP
	for j = 1:length(alldata)
		y = mean(alldata(j).LFP_xcorr(:,:,i),2);
		y = y/max(y);
		plot(-chunk_size+1:chunk_size-1,y,'Color',c(j,:))
	end
end

subplot(2,3,4); hold on
for i = 1:length(alldata)
	plot(mean(alldata(i).PID(a:z,:)),alldata(i).LFP_lags,'+-','Color',c(i,:))
end
xlabel('\mu_{Stimulus}')
ylabel('Stim \rightarrow LFP Lag (ms)')
legend({'Low Ca','Normal Ca','High Ca'})

prettyFig()

if being_published
	snapnow
	delete(gcf)
end

%% Version Info
%
pFooter;



