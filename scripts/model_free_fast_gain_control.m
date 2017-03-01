

%% Fast gain control in naturalistic stimulus
% In this document, I attempt to use a model-free analysis method to demonstrate that there exists fast gain control in the naturalistic stimulus. 



pHeader;

%% Analysis pipeline 
% This method of analysis is so fiendishly complicated that I made a handy graphic that explains what I'm doing: 

figure('outerposition',[0 0 1000 655],'PaperUnits','points','PaperSize',[1000 655]); hold on
% show the cartoons
o = imread('../images/fast_gain_control_analysis.png');
imagesc(o);
axis ij
axis image
axis off

prettyFig();

if being_published
	snapnow
	delete(gcf)
end



cdata = consolidateData2(getPath(dataManager,'4608c42b12191b383c84fea52392ea97'));
[cdata, data] =  assembleScaledNatStim(cdata);


;;;;;;;;  ;;;;;;;;    ;;;    ;;          ;;;;;;;;     ;;;    ;;;;;;;;    ;;;    
;;     ;; ;;         ;; ;;   ;;          ;;     ;;   ;; ;;      ;;      ;; ;;   
;;     ;; ;;        ;;   ;;  ;;          ;;     ;;  ;;   ;;     ;;     ;;   ;;  
;;;;;;;;  ;;;;;;   ;;     ;; ;;          ;;     ;; ;;     ;;    ;;    ;;     ;; 
;;   ;;   ;;       ;;;;;;;;; ;;          ;;     ;; ;;;;;;;;;    ;;    ;;;;;;;;; 
;;    ;;  ;;       ;;     ;; ;;          ;;     ;; ;;     ;;    ;;    ;;     ;; 
;;     ;; ;;;;;;;; ;;     ;; ;;;;;;;;    ;;;;;;;;  ;;     ;;    ;;    ;;     ;; 

%% Real data: ab2A responses to 2-butanone
% In this section I apply this analysis method to ab2A responses to naturalistic stimuli constructed with 2-butanone. Since this stimulus is presented many times in scaled versions of itself, it covers the entire dynamic range of the neuron. (a-c) Analysis of fast gain control in LFP. (a) LFP responses vs. whiff amplitude, for sets of whiffs with similar amplitude (see methods in previous figure). The variation in whiff amplitude cannot describe the variation in the LFP responses. Lines connect whiffs of the same size that come from the same set. (b) LFP responses vs. mean stimulus in the preceding 300ms. The strong negative correlation shows what we saw in the raw data: that LFP responses to whiffs following other whiffs are smaller than LFp responses to whiffs that occur in isolation. (c) Spearman correlation between LFP responses and whiff amplitude (black) and between LFP responses and mean stimulus in the preceding window. Crosses indicate significance (p<0.05). Note that there is a pronounced dip at around ~300ms. 

%%
% (d-f) Similar analysis, but for the firing rates. 

[s,x,r,whiff_s,group] = binnedWhiffGainAnalysis(data(2).S,data(2).X,data(2).R);

figure('outerposition',[0 0 1400 901],'PaperUnits','points','PaperSize',[1400 901]); hold on
subplot(2,3,1); hold on
for i = 1:max(group)
	tempx = whiff_s(group==i)/mean(whiff_s(group==i));
	tempy = x(group==i);
	[tempx,idx] = sort(tempx);
	plot(tempx,tempy(idx),'+-');
end
l = plot(NaN,NaN,'k+');
legend(l,['\rho = ' oval(spear(whiff_s,x))],'Location','northeast')
xlabel('Whiff amplitude (norm)')
ylabel('LFP response (norm)')

subplot(2,3,2); hold on
for i = 1:max(group)
	tempx = s(group==i);
	tempy = x(group==i);
	[tempx,idx] = sort(tempx);
	plot(tempx,tempy(idx),'+-');
end
l = plot(NaN,NaN,'k+');
legend(l,['\rho = ' oval(spear(s,x))],'Location','northeast')
xlabel('\mu_{Stimulus} in preceding 300ms')
ylabel('LFP response (norm)')

subplot(2,3,4); hold on
for i = 1:max(group)
	tempx = whiff_s(group==i)/mean(whiff_s(group==i));
	tempy = r(group==i);
	[tempx,idx] = sort(tempx);
	plot(tempx,tempy(idx),'+-');
end
l = plot(NaN,NaN,'k+');
legend(l,['\rho = ' oval(spear(whiff_s,r))],'Location','northeast')
xlabel('Whiff amplitude (norm)')
ylabel('Firing response (norm)')

subplot(2,3,5); hold on
for i = 1:max(group)
	tempx = s(group==i);
	tempy = r(group==i);
	[tempx,idx] = sort(tempx);
	plot(tempx,tempy(idx),'+-');
end
l = plot(NaN,NaN,'k+');
legend(l,['\rho = ' oval(spear(s,r))],'Location','northeast')
xlabel('\mu_{Stimulus} in preceding 300ms')
ylabel('Firing response (norm)')

% now vary the stimulus history length 
all_history_lengths = unique(floor(logspace(1,4,60)));

rho_x_whiff = NaN*all_history_lengths;
rho_x_gain = NaN*all_history_lengths;
rho_r_whiff = NaN*all_history_lengths;
rho_r_gain = NaN*all_history_lengths;

rho_x_whiff_p = NaN*all_history_lengths;
rho_x_gain_p = NaN*all_history_lengths;
rho_r_whiff_p = NaN*all_history_lengths;
rho_r_gain_p = NaN*all_history_lengths;

for i = 1:length(all_history_lengths)
	[s,x,r,whiff_s,group] = binnedWhiffGainAnalysis(data(2).S,data(2).X,data(2).R,'history_length',all_history_lengths(i));

	% compute correlations 
	[rho_x_gain(i), rho_x_gain_p(i)] = corr(s,x,'type','Spearman');
	[rho_x_whiff(i), rho_x_whiff_p(i)] = corr(whiff_s,x,'type','Spearman');
	[rho_r_gain(i), rho_r_gain_p(i)] = corr(s,r,'type','Spearman');
	[rho_r_whiff(i), rho_r_whiff_p(i)] = corr(whiff_s,r,'type','Spearman');


end

% emphasise significant points
subplot(2,3,3); hold on
plot(all_history_lengths,rho_x_whiff,'k');
plot(all_history_lengths,rho_x_gain,'r');


subplot(2,3,6); hold on
plot(all_history_lengths,rho_r_whiff,'k');
plot(all_history_lengths,rho_r_gain,'r');



rho_x_gain(rho_x_gain_p>.05) = NaN;
rho_x_whiff(rho_x_whiff_p>.05) = NaN;
rho_r_gain(rho_r_gain_p>.05) = NaN;
rho_r_whiff(rho_r_whiff_p>.05) = NaN;

subplot(2,3,3); hold on
plot(all_history_lengths,rho_x_whiff,'k+');
plot(all_history_lengths,rho_x_gain,'r+');
set(gca,'XScale','log')
xlabel('Stimulus history length (ms)')
ylabel('\rho')

subplot(2,3,6); hold on
plot(all_history_lengths,rho_r_whiff,'k+');
plot(all_history_lengths,rho_r_gain,'r+');
set(gca,'XScale','log')
xlabel('Stimulus history length (ms)')
ylabel('\rho')


prettyFig();
labelFigure('x_offset',0)

if being_published
	snapnow
	delete(gcf)
end

%%
% Now I make another plot where I plot 

;;     ;;    ;;;    ;;       ;;;; ;;;;;;;;     ;;;    ;;;;;;;; ;;;;  ;;;;;;;  ;;    ;; 
;;     ;;   ;; ;;   ;;        ;;  ;;     ;;   ;; ;;      ;;     ;;  ;;     ;; ;;;   ;; 
;;     ;;  ;;   ;;  ;;        ;;  ;;     ;;  ;;   ;;     ;;     ;;  ;;     ;; ;;;;  ;; 
;;     ;; ;;     ;; ;;        ;;  ;;     ;; ;;     ;;    ;;     ;;  ;;     ;; ;; ;; ;; 
 ;;   ;;  ;;;;;;;;; ;;        ;;  ;;     ;; ;;;;;;;;;    ;;     ;;  ;;     ;; ;;  ;;;; 
  ;; ;;   ;;     ;; ;;        ;;  ;;     ;; ;;     ;;    ;;     ;;  ;;     ;; ;;   ;;; 
   ;;;    ;;     ;; ;;;;;;;; ;;;; ;;;;;;;;  ;;     ;;    ;;    ;;;;  ;;;;;;;  ;;    ;; 

;;    ;; ;;          ;;     ;;  ;;;;;;;  ;;;;;;;;  ;;;;;;;; ;;       
;;;   ;; ;;          ;;;   ;;; ;;     ;; ;;     ;; ;;       ;;       
;;;;  ;; ;;          ;;;; ;;;; ;;     ;; ;;     ;; ;;       ;;       
;; ;; ;; ;;          ;; ;;; ;; ;;     ;; ;;     ;; ;;;;;;   ;;       
;;  ;;;; ;;          ;;     ;; ;;     ;; ;;     ;; ;;       ;;       
;;   ;;; ;;          ;;     ;; ;;     ;; ;;     ;; ;;       ;;       
;;    ;; ;;;;;;;;    ;;     ;;  ;;;;;;;  ;;;;;;;;  ;;;;;;;; ;;;;;;;; 


%% Validation: NL model (LFP)
% As a validation, I repeat this analysis using synthetic data generated by a NL model. (a-c) Analysis using synthetic data generated by a NL model fit to the LFP responses. (d-f) Analysis using real LFP responses. (a) NL model responses vs. whiff amplitude, for similar sized whiffs. Whiff amplitudes cannot account for observed variation in responses. (b) NL model responses vs. mean stimulus in preceding 300ms. The mean stimulus in the preceding 300ms can account for the variation in the synthetic data, but it goes the "wrong way" -- NL model responses tend to be bigger when the mean stimulus in the preceding 300ms is larger -- consistent with an integrating, but non-adapting system. (c) Spearman correlation for NL model responses vs. whiff amplitude (black) and vs. mean stimulus in preceding window (red). Crossed indicate significant correlation. Note that for this synthetic data set, the mean stimulus in any preceding window always tends to amplify responses. 

%%
% (e-f) as in the previous figure. 

clear p
p.k_D = 0.21246;
p.n = 1.4688;


% generate responses using this model 
warning off
i = 2;
for j = 1:size(data(i).X,2)
	data(i).P(:,j) = NLModel([data(i).S(:,j) - min(data(i).S(:,j)) data(i).X(:,j)] ,p);
end

warning on

% we now attempt to estimate a gain control timescale from the data
% first collect all whiffs
[s,x,r,whiff_s,group] = binnedWhiffGainAnalysis(data(2).S,data(2).P,-data(2).X);

figure('outerposition',[0 0 1400 901],'PaperUnits','points','PaperSize',[1400 901]); hold on
subplot(2,3,1); hold on
for i = 1:max(group)
	tempx = whiff_s(group==i)/mean(whiff_s(group==i));
	tempy = x(group==i);
	[tempx,idx] = sort(tempx);
	plot(tempx,tempy(idx),'+-');
end
l = plot(NaN,NaN,'k+');
legend(l,['\rho = ' oval(spear(whiff_s,x))],'Location','northeast')
xlabel('Whiff amplitude (norm)')
ylabel('NL model response (norm)')

subplot(2,3,2); hold on
for i = 1:max(group)
	tempx = s(group==i);
	tempy = x(group==i);
	[tempx,idx] = sort(tempx);
	plot(tempx,tempy(idx),'+-');
end
l = plot(NaN,NaN,'k+');
legend(l,['\rho = ' oval(spear(s,x))],'Location','northeast')
xlabel('\mu_{Stimulus} in preceding 300ms')
ylabel('NL model response (norm)')

subplot(2,3,4); hold on
for i = 1:max(group)
	tempx = whiff_s(group==i)/mean(whiff_s(group==i));
	tempy = r(group==i);
	[tempx,idx] = sort(tempx);
	plot(tempx,tempy(idx),'+-');
end
l = plot(NaN,NaN,'k+');
legend(l,['\rho = ' oval(spear(whiff_s,r))],'Location','northeast')
xlabel('Whiff amplitude (norm)')
ylabel('LFP response (norm)')

subplot(2,3,5); hold on
for i = 1:max(group)
	tempx = s(group==i);
	tempy = r(group==i);
	[tempx,idx] = sort(tempx);
	plot(tempx,tempy(idx),'+-');
end
l = plot(NaN,NaN,'k+');
legend(l,['\rho = ' oval(spear(s,r))],'Location','northeast')
xlabel('\mu_{Stimulus} in preceding 300ms')
ylabel('LFP response (norm)')

% now vary the stimulus history length 
all_history_lengths = unique(floor(logspace(1,4,60)));

rho_x_whiff = NaN*all_history_lengths;
rho_x_gain = NaN*all_history_lengths;
rho_r_whiff = NaN*all_history_lengths;
rho_r_gain = NaN*all_history_lengths;

rho_x_whiff_p = NaN*all_history_lengths;
rho_x_gain_p = NaN*all_history_lengths;
rho_r_whiff_p = NaN*all_history_lengths;
rho_r_gain_p = NaN*all_history_lengths;

for i = 1:length(all_history_lengths)
	[s,x,r,whiff_s,group] = binnedWhiffGainAnalysis(data(2).S,data(2).P,-data(2).X,'history_length',all_history_lengths(i));

	% compute correlations 
	[rho_x_gain(i), rho_x_gain_p(i)] = corr(s,x,'type','Spearman');
	[rho_x_whiff(i), rho_x_whiff_p(i)] = corr(whiff_s,x,'type','Spearman');
	[rho_r_gain(i), rho_r_gain_p(i)] = corr(s,r,'type','Spearman');
	[rho_r_whiff(i), rho_r_whiff_p(i)] = corr(whiff_s,r,'type','Spearman');


end

% emphasise significant points
subplot(2,3,3); hold on
plot(all_history_lengths,rho_x_whiff,'k');
plot(all_history_lengths,rho_x_gain,'r');


subplot(2,3,6); hold on
plot(all_history_lengths,rho_r_whiff,'k');
plot(all_history_lengths,rho_r_gain,'r');



rho_x_gain(rho_x_gain_p>.05) = NaN;
rho_x_whiff(rho_x_whiff_p>.05) = NaN;
rho_r_gain(rho_r_gain_p>.05) = NaN;
rho_r_whiff(rho_r_whiff_p>.05) = NaN;

subplot(2,3,3); hold on
plot(all_history_lengths,rho_x_whiff,'k+');
plot(all_history_lengths,rho_x_gain,'r+');
set(gca,'XScale','log')
xlabel('Stimulus history length (ms)')
ylabel('\rho')

subplot(2,3,6); hold on
plot(all_history_lengths,rho_r_whiff,'k+');
plot(all_history_lengths,rho_r_gain,'r+');
set(gca,'XScale','log')
xlabel('Stimulus history length (ms)')
ylabel('\rho')


prettyFig();
labelFigure('x_offset',0);

if being_published
	snapnow
	delete(gcf)
end


;;     ;;    ;;;    ;;       ;;;; ;;;;;;;;     ;;;    ;;;;;;;; ;;;;  ;;;;;;;  ;;    ;; 
;;     ;;   ;; ;;   ;;        ;;  ;;     ;;   ;; ;;      ;;     ;;  ;;     ;; ;;;   ;; 
;;     ;;  ;;   ;;  ;;        ;;  ;;     ;;  ;;   ;;     ;;     ;;  ;;     ;; ;;;;  ;; 
;;     ;; ;;     ;; ;;        ;;  ;;     ;; ;;     ;;    ;;     ;;  ;;     ;; ;; ;; ;; 
 ;;   ;;  ;;;;;;;;; ;;        ;;  ;;     ;; ;;;;;;;;;    ;;     ;;  ;;     ;; ;;  ;;;; 
  ;; ;;   ;;     ;; ;;        ;;  ;;     ;; ;;     ;;    ;;     ;;  ;;     ;; ;;   ;;; 
   ;;;    ;;     ;; ;;;;;;;; ;;;; ;;;;;;;;  ;;     ;;    ;;    ;;;;  ;;;;;;;  ;;    ;; 

;;    ;; ;;       ;;    ;;    ;;     ;;  ;;;;;;;  ;;;;;;;;  ;;;;;;;; ;;       
;;;   ;; ;;       ;;;   ;;    ;;;   ;;; ;;     ;; ;;     ;; ;;       ;;       
;;;;  ;; ;;       ;;;;  ;;    ;;;; ;;;; ;;     ;; ;;     ;; ;;       ;;       
;; ;; ;; ;;       ;; ;; ;;    ;; ;;; ;; ;;     ;; ;;     ;; ;;;;;;   ;;       
;;  ;;;; ;;       ;;  ;;;;    ;;     ;; ;;     ;; ;;     ;; ;;       ;;       
;;   ;;; ;;       ;;   ;;;    ;;     ;; ;;     ;; ;;     ;; ;;       ;;       
;;    ;; ;;;;;;;; ;;    ;;    ;;     ;;  ;;;;;;;  ;;;;;;;;  ;;;;;;;; ;;;;;;;; 


%% Validation: NLN model (firing rate)
% As a validation, I repeat this analysis using synthetic data generated by a NLN model. This NLN model has a negative lobe, since it is fit to the firing rate data, and this could in principle give us something that looks like "fast gain control".  (a-c) Analysis using synthetic data generated by a NLN model fit to the firing rate responses. (d-f) Analysis using real firing rate responses. (a) NLN model responses vs. whiff amplitude, for similar sized whiffs. Whiff amplitudes cannot account for observed variation in responses. (b) NLN model responses vs. mean stimulus in preceding 300ms. The mean stimulus in the preceding 300ms can account for the variation in the synthetic data, similar to what we see in the data. (c) Spearman correlation for NLN model responses vs. whiff amplitude (black) and vs. mean stimulus in preceding window (red). Crossed indicate significant correlation. Note that for this synthetic data set, the mean stimulus in any preceding window always tends to amplify responses. 

%%
% (d-f) as in the previous figure with real data. 

clear p
p.  k_D = 0.1109;
p.    n = .9688;

% generate responses using this model 
warning off
i = 2;
for j = 1:size(data(i).X,2)
	try
		data(i).P(:,j) = NLNmodel([data(i).S(:,j) - min(data(i).S(:,j)) data(i).R(:,j)] ,p);
	catch
	end
end

warning on

% we now attempt to estimate a gain control timescale from the data
% first collect all whiffs
[s,x,r,whiff_s,group] = binnedWhiffGainAnalysis(data(2).S,-data(2).P,data(2).R);

figure('outerposition',[0 0 1400 901],'PaperUnits','points','PaperSize',[1400 901]); hold on
subplot(2,3,1); hold on
for i = 1:max(group)
	tempx = whiff_s(group==i)/mean(whiff_s(group==i));
	tempy = x(group==i);
	[tempx,idx] = sort(tempx);
	plot(tempx,tempy(idx),'+-');
end
l = plot(NaN,NaN,'k+');
legend(l,['\rho = ' oval(spear(whiff_s,x))],'Location','northeast')
xlabel('Whiff amplitude (norm)')
ylabel('NLN model response (norm)')

subplot(2,3,2); hold on
for i = 1:max(group)
	tempx = s(group==i);
	tempy = x(group==i);
	[tempx,idx] = sort(tempx);
	plot(tempx,tempy(idx),'+-');
end
l = plot(NaN,NaN,'k+');
legend(l,['\rho = ' oval(spear(s,x))],'Location','northeast')
xlabel('\mu_{Stimulus} in preceding 300ms')
ylabel('NLN model response (norm)')

subplot(2,3,4); hold on
for i = 1:max(group)
	tempx = whiff_s(group==i)/mean(whiff_s(group==i));
	tempy = r(group==i);
	[tempx,idx] = sort(tempx);
	plot(tempx,tempy(idx),'+-');
end
l = plot(NaN,NaN,'k+');
legend(l,['\rho = ' oval(spear(whiff_s,r))],'Location','northeast')
xlabel('Whiff amplitude (norm)')
ylabel('Firing response (norm)')

subplot(2,3,5); hold on
for i = 1:max(group)
	tempx = s(group==i);
	tempy = r(group==i);
	[tempx,idx] = sort(tempx);
	plot(tempx,tempy(idx),'+-');
end
l = plot(NaN,NaN,'k+');
legend(l,['\rho = ' oval(spear(s,r))],'Location','northeast')
xlabel('\mu_{Stimulus} in preceding 300ms')
ylabel('Firing response (norm)')

% now vary the stimulus history length 
all_history_lengths = unique(floor(logspace(1,4,60)));

rho_x_whiff = NaN*all_history_lengths;
rho_x_gain = NaN*all_history_lengths;
rho_r_whiff = NaN*all_history_lengths;
rho_r_gain = NaN*all_history_lengths;

rho_x_whiff_p = NaN*all_history_lengths;
rho_x_gain_p = NaN*all_history_lengths;
rho_r_whiff_p = NaN*all_history_lengths;
rho_r_gain_p = NaN*all_history_lengths;

for i = 1:length(all_history_lengths)
	[s,x,r,whiff_s,group] = binnedWhiffGainAnalysis(data(2).S,-data(2).P,data(2).R,'history_length',all_history_lengths(i));

	% compute correlations 
	[rho_x_gain(i), rho_x_gain_p(i)] = corr(s,x,'type','Spearman');
	[rho_x_whiff(i), rho_x_whiff_p(i)] = corr(whiff_s,x,'type','Spearman');
	[rho_r_gain(i), rho_r_gain_p(i)] = corr(s,r,'type','Spearman');
	[rho_r_whiff(i), rho_r_whiff_p(i)] = corr(whiff_s,r,'type','Spearman');


end

% emphasise significant points
subplot(2,3,3); hold on
plot(all_history_lengths,rho_x_whiff,'k');
plot(all_history_lengths,rho_x_gain,'r');


subplot(2,3,6); hold on
plot(all_history_lengths,rho_r_whiff,'k');
plot(all_history_lengths,rho_r_gain,'r');



rho_x_gain(rho_x_gain_p>.05) = NaN;
rho_x_whiff(rho_x_whiff_p>.05) = NaN;
rho_r_gain(rho_r_gain_p>.05) = NaN;
rho_r_whiff(rho_r_whiff_p>.05) = NaN;

subplot(2,3,3); hold on
plot(all_history_lengths,rho_x_whiff,'k+');
plot(all_history_lengths,rho_x_gain,'r+');
set(gca,'XScale','log')
xlabel('Stimulus history length (ms)')
ylabel('\rho')

subplot(2,3,6); hold on
plot(all_history_lengths,rho_r_whiff,'k+');
plot(all_history_lengths,rho_r_gain,'r+');
set(gca,'XScale','log')
xlabel('Stimulus history length (ms)')
ylabel('\rho')


prettyFig();

labelFigure('x_offset',0)

if being_published
	snapnow
	delete(gcf)
end

%% Version Info
%
pFooter;


