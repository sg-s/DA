%%
%  In this document, I fit a single adapting NLN model (aNLN4) to all the data.

pHeader;


% get the filter from the Gaussian stimuli 
clear MSGdata
MSGdata = consolidateData2(getPath(dataManager,'93ba5d68174e3df9f462a1fc48c581da'));
MSGdata = cleanMSGdata(MSGdata);

% make sure stimulus is always positive
for i = 1:size(MSGdata.PID,2)
	MSGdata.PID(:,i) = MSGdata.PID(:,i) - min(MSGdata.PID(:,i));
end

% a = 25e3; z = 55e3;
% clear data
% for i = 1:max(MSGdata.paradigm)
% 	S = MSGdata.PID(:,MSGdata.paradigm==i);
% 	R = MSGdata.fA(:,MSGdata.paradigm==i);
% 	rm_this = sum(R) == 0;
% 	S(:,rm_this) = []; R(:,rm_this) = [];
% 	data(i).stimulus = mean(S(a:z,:),2);
% 	data(i).response = mean(R(a:z,:),2);
% 	data(i).response(1:10e3) = NaN;
% end

;;     ;;  ;;;;;;;  ;;;;;;;;  ;;;;;;;; ;;       
;;;   ;;; ;;     ;; ;;     ;; ;;       ;;       
;;;; ;;;; ;;     ;; ;;     ;; ;;       ;;       
;; ;;; ;; ;;     ;; ;;     ;; ;;;;;;   ;;       
;;     ;; ;;     ;; ;;     ;; ;;       ;;       
;;     ;; ;;     ;; ;;     ;; ;;       ;;       
;;     ;;  ;;;;;;;  ;;;;;;;;  ;;;;;;;; ;;;;;;;; 

;;;;;;;;     ;;;    ;;;;;;;;     ;;;    ;;     ;;  ;;;;;;  
;;     ;;   ;; ;;   ;;     ;;   ;; ;;   ;;;   ;;; ;;    ;; 
;;     ;;  ;;   ;;  ;;     ;;  ;;   ;;  ;;;; ;;;; ;;       
;;;;;;;;  ;;     ;; ;;;;;;;;  ;;     ;; ;; ;;; ;;  ;;;;;;  
;;        ;;;;;;;;; ;;   ;;   ;;;;;;;;; ;;     ;;       ;; 
;;        ;;     ;; ;;    ;;  ;;     ;; ;;     ;; ;;    ;; 
;;        ;;     ;; ;;     ;; ;;     ;; ;;     ;;  ;;;;;;  


clear p
p.    k0 = 0.0968;
p.     B = 0.8289;
p.tau_y1 = 27.1055;
p.tau_y2 = 96.9453;
p.     A = 0.7000;
p.     C = 178.9440;
p.     n = 2;


figure('outerposition',[0 0 1101 901],'PaperUnits','points','PaperSize',[1101 901]); hold on
clear ax
ax(1) = subplot(4,3,1); hold on
ax(2) = subplot(4,3,2); hold on
ax(3) = subplot(4,3,3); hold on


% first, show 3 nonlinearities from the Gaussians with the right colours 
c = parula(11);
for i = [2 6 10]
	this_paradigm = find(MSGdata.paradigm == i,1,'first');
	S = MSGdata.PID(:,this_paradigm);
	[R, a, b, kD, K] = aNLN4(S,p);
	% plot actual curve with mean kD
	x = logspace(-2,2,100);
	A = 1./(1 + mean(kD(35e3:55e3))./x);
	plot(ax(1),x,A,'Color',c(i,:))
	% now plot the points on it that are actually achieved 
	plot(ax(1),S(35e3:55e3),a(35e3:55e3),'.','Color',c(i,:))
end
set(ax(1),'XScale','log','XLim',[1e-2 10])
ax(1).YAxisLocation = 'right';
xlabel(ax(1),'Stimulus (V)')
ylabel(ax(1),'a')
axes(ax(1))
axis square

% show the filter 
plot(ax(2),K/norm(K),'k')
set(ax(2),'XLim',[0 600])
xlabel(ax(2),'Lag (ms)')
ylabel(ax(2),'Filter (norm)')
axes(ax(2))
axis square

;;     ;;  ;;;;;;   ;;;;;;   
;;;   ;;; ;;    ;; ;;    ;;  
;;;; ;;;; ;;       ;;        
;; ;;; ;;  ;;;;;;  ;;   ;;;; 
;;     ;;       ;; ;;    ;;  
;;     ;; ;;    ;; ;;    ;;  
;;     ;;  ;;;;;;   ;;;;;;   

a = 35e3;
z = 55e3;

% generate responses using the model 
if exist('.cache/aNLN4_responses_to_MSG.mat','file') == 0
	NLN_pred = NaN*MSGdata.PID;
	for i = 1:length(MSGdata.paradigm)
		S = MSGdata.PID(:,i) - min(MSGdata.PID(:,i));
		NLN_pred(:,i) = aNLN4(S,p);
	end
	save('.cache/aNLN4_responses_to_MSG.mat','NLN_pred')
else
	load('.cache/aNLN4_responses_to_MSG.mat')
end
MSGdata.NLN_pred = NLN_pred;
clear NLN_pred

% back out new linear filters for this
time = 1e-3*(1:length(MSGdata.PID));
if exist('.cache/aNLN4_MSG_linear_prediction.mat','file') == 0
	NLN_fp = NaN*MSGdata.PID;
	for i = 1:length(MSGdata.paradigm)
		S = MSGdata.PID(:,i);
		R = MSGdata.NLN_pred(:,i);
		K = fitFilter2Data(S(a:z),R(a:z),'filter_length',1e3,'offset',200);
		K = K(100:end-100);
		filtertime = 1e-3*(1:length(K)) - .1;
		NLN_fp(:,i) = convolve(time,S,K,filtertime);
		K = K/(nanstd(NLN_fp(a:z,i))/nanstd(S(a:z))); % normalise correctly 
		NLN_fp(:,i) = convolve(time,S,K,filtertime);
	end
	save('.cache/aNLN4_MSG_linear_prediction.mat','NLN_fp')
else
	load('.cache/aNLN4_MSG_linear_prediction.mat')
end
MSGdata.NLN_fp = NLN_fp;
clear NLN_fp

% recreate the I/O plot 
subplot(4,3,7); hold on; cla
ss = 100;
all_x = 0:0.1:2;
a = 35e3; z = 55e3;
c = parula(10);
for i = 1:max(MSGdata.paradigm) % iterate over all paradigms 
	y = MSGdata.NLN_pred(a:z,MSGdata.paradigm == i);
	x = MSGdata.NLN_fp(a:z,MSGdata.paradigm == i);
	s = MSGdata.PID(a:z,MSGdata.paradigm == i);
	rm_this = sum(y)==0;
	y(:,rm_this) = [];
	x(:,rm_this) = [];
	s(:,rm_this) = [];
	y = nanmean(y,2);
	x = nanmean(x,2);
	s = nanmean(s,2);
	x = x - nanmean(x);
	x = x + nanmean(nanmean(s));
	plotPieceWiseLinear(x,y,'nbins',50,'Color',c(i,:),'show_error',false,'LineWidth',3);
end
xlabel('Projected stimulus (V)')
ylabel('NLN response (Hz)')

% compute gains and plot them
MSGdata.NLN_gain = NaN*MSGdata.paradigm;
for i = 1:length(MSGdata.paradigm)
	x = MSGdata.NLN_fp(a:z,i);
	y = MSGdata.NLN_pred(a:z,i);
	try
		ff = fit(x(:),y(:),'poly1');
		MSGdata.NLN_gain(i) = ff.p1;
	catch
	end
end

subplot(4,3,8); hold on; cla
mean_stim = nanmean(MSGdata.PID(a:z,:));

% correct for some trivial scaling 
x = MSGdata.NLN_gain(:);
y = MSGdata.fA_gain(:);
rm_this = isnan(x) | isnan(y);
ff = fit(x(~rm_this),y(~rm_this),'poly1');
MSGdata.NLN_gain = ff(MSGdata.NLN_gain);

% show gain changes -- gain vs. mean stimulus
all_x = []; all_y = [];
for i = 1:max(MSGdata.paradigm)
	x = mean_stim(MSGdata.paradigm==i);
	y = MSGdata.NLN_gain(MSGdata.paradigm==i);
	rm_this = y == 0 | isnan(y);
	x(rm_this) = []; y(rm_this) = []; 
	all_x = [all_x(:); x(:)]; all_y = [all_y(:); y(:)];
	plot(x,y,'+','Color',c(i,:));
end

% fit a power law with exponent -1
mean_stim = mean_stim(:);
g = MSGdata.NLN_gain(:);
options = fitoptions(fittype('power1'));
options.Lower = [-Inf -1];
options.Upper = [Inf -1];
cf = fit(mean_stim(~isnan(g)),g(~isnan(g)),'power1',options);
plot(sort(mean_stim),cf(sort(mean_stim)),'r');
set(gca,'XScale','log','YScale','log','YLim',[10 300],'XLim',[.1 2.5])
xlabel('\mu_{Stimulus} (V)')
ylabel('NLN gain (Hz/V)')

% compare model gains to actual gains
subplot(4,3,9); hold on; cla
for i = 1:max(MSGdata.paradigm)
	plot(MSGdata.NLN_gain(MSGdata.paradigm == i),MSGdata.fA_gain(MSGdata.paradigm == i),'+','Color',c(i,:))
end
clear l
l = plot(NaN,NaN,'k+');
legend(l,['r^2 = ' oval(rsquare(MSGdata.NLN_gain,MSGdata.fA_gain))],'Location','southeast')
xlabel('NLN gain (Hz/V)')
ylabel('Observed gain (Hz/V)')
plot([0 250],[0 250],'k--')


;;    ;;    ;;;    ;;;;;;;; ;;     ;; ;;;;;;;;     ;;;    ;;       
;;;   ;;   ;; ;;      ;;    ;;     ;; ;;     ;;   ;; ;;   ;;       
;;;;  ;;  ;;   ;;     ;;    ;;     ;; ;;     ;;  ;;   ;;  ;;       
;; ;; ;; ;;     ;;    ;;    ;;     ;; ;;;;;;;;  ;;     ;; ;;       
;;  ;;;; ;;;;;;;;;    ;;    ;;     ;; ;;   ;;   ;;;;;;;;; ;;       
;;   ;;; ;;     ;;    ;;    ;;     ;; ;;    ;;  ;;     ;; ;;       
;;    ;; ;;     ;;    ;;     ;;;;;;;  ;;     ;; ;;     ;; ;;;;;;;; 

 ;;;;;;  ;;;;;;;; ;;;; ;;     ;; 
;;    ;;    ;;     ;;  ;;;   ;;; 
;;          ;;     ;;  ;;;; ;;;; 
 ;;;;;;     ;;     ;;  ;; ;;; ;; 
      ;;    ;;     ;;  ;;     ;; 
;;    ;;    ;;     ;;  ;;     ;; 
 ;;;;;;     ;;    ;;;; ;;     ;; 


% show the naturalistic stimulus fits

cdata = consolidateData2(getPath(dataManager,'4608c42b12191b383c84fea52392ea97'));
[cdata, data] =  assembleScaledNatStim(cdata);
time = 1e-3*(1:length(data(1).S));

% show the time series of the natualistic stimulus
ax(4) = subplot(4,3,4:5); hold on
clear l
l(1) = plot(ax(4),time,data(2).R(:,3),'k');
% generate responses and plot them too

for j = 1:size(data(2).X,2)
	data(2).P(:,j) = aNLN4(data(2).S(:,j) - min(data(2).S(:,j)) ,p);
end
% fix some trivial scaling
ff = fit(data(2).P(:),data(2).R(:),'poly1');
for j = 1:size(data(2).X,2)
	data(2).P(:,j) = ff(data(2).P(:,j));
end

l(2) = plot(ax(4),time,data(2).P(:,3),'r');
legend(l,{'ab2A',['model r^2 = ' oval(rsquare(data(2).P(:,3),data(2).R(:,3)))]})

xlabel(ax(4),'Time (s)')
ylabel(ax(4),'Firing rate (Hz)')

% show r^2 for every whiff
ax(5) = subplot(4,3,6); hold on
x = []; y = [];
for i = 2
	for j = 1:3
		S = data(i).S(:,j);
		P = data(i).P(:,j);
		R = data(i).R(:,j);
		ws = whiffStatistics(S,R,R,300,'MinPeakProminence',max(S/1e2),'debug',false);
		y = [y; ws.peak_firing_rate];
		ws = whiffStatistics(S,P,P,300,'MinPeakProminence',max(S/1e2),'debug',false);
		x = [x; ws.peak_firing_rate];
	end
end
clear l
l = plot(ax(5),x,y,'k.','MarkerSize',20);
legend(l,['r^2 = ' oval(rsquare(x,y))],'Location','southeast')
xlabel(ax(5),'NLN model response (Hz)')
ylabel(ax(5),'ab2A response (Hz)')
title(ax(5),'Whiff-specific responses')
set(ax(5),'XLim',[0 300],'YLim',[0 300])


prettyFig('fs',12);

if being_published
	snapnow
	delete(gcf)
end

return


;;     ;;    ;;;    ;;;;;;;;  ;;;;    ;;;    ;;    ;;  ;;;;;;  ;;;;;;;; 
;;     ;;   ;; ;;   ;;     ;;  ;;    ;; ;;   ;;;   ;; ;;    ;; ;;       
;;     ;;  ;;   ;;  ;;     ;;  ;;   ;;   ;;  ;;;;  ;; ;;       ;;       
;;     ;; ;;     ;; ;;;;;;;;   ;;  ;;     ;; ;; ;; ;; ;;       ;;;;;;   
 ;;   ;;  ;;;;;;;;; ;;   ;;    ;;  ;;;;;;;;; ;;  ;;;; ;;       ;;       
  ;; ;;   ;;     ;; ;;    ;;   ;;  ;;     ;; ;;   ;;; ;;    ;; ;;       
   ;;;    ;;     ;; ;;     ;; ;;;; ;;     ;; ;;    ;;  ;;;;;;  ;;;;;;;; 


% get the variance switching data
clear VSdata
[VSdata.PID, VSdata.LFP, VSdata.fA, VSdata.paradigm, VSdata.orn, VSdata.fly] = consolidateData(getPath(dataManager,'e30707e8e8ef6c0d832eee31eaa585aa'),1);
% remove baseline from stimulus
VSdata.PID = bsxfun(@minus, VSdata.PID, min(VSdata.PID));
VSdata.LFP = bsxfun(@minus, VSdata.LFP, mean(VSdata.LFP(1:4000,:)));

% get the variance switching data filters 
load(getPath(dataManager,'457ee16a326f47992e35a7d5281f9cc4'));
VSdata.K1 = nanmean(K1,2);


% % this fits the data
% this_trial = 2;
% clear data
% data.response = -LFP(50e3:90e3,this_trial);
% data.stimulus = PID(50e3:90e3,this_trial);
% data.response(1:10e3) = NaN;

PID = VSdata.PID;
LFP = VSdata.LFP;

global_start = 40e3; % 40 seconds
global_end = length(PID) - 15e3; 

% generate responses using the LFP model
clear p
p.  n  =2.1805;
p.   A = 0.5680;
p.tauA = 31.7382;
p.tauB = 124.2783;
p.   C = 0.9932;
p.   D = -9.7793;
p.  k0 = 0.3723;
p.   B = 0.1006;

NL_pred = NaN*LFP;
for i = 1:width(LFP)
	NL_pred(:,i) = aNL4(PID(:,i),p);
end



% filter to remove spikes
for i = 1:width(LFP)
	LFP(:,i) = filtfilt(ones(30,1),30,LFP(:,i));
end


% reshape the LFP signals
block_length = 1e4;
reshaped_LFP = LFP(global_start:end-1e4-1,1:width(PID));
reshaped_LFP = reshape(reshaped_LFP,block_length,width(reshaped_LFP)*length(reshaped_LFP)/block_length);

% also reshape the PID
reshaped_PID = PID(global_start:end-1e4-1,1:width(PID));
reshaped_PID = reshape(reshaped_PID,block_length,width(reshaped_PID)*length(reshaped_PID)/block_length);

% also reshape the prediction
reshaped_NL_pred = NL_pred(global_start:end-1e4-1,1:width(NL_pred));
reshaped_NL_pred = reshape(reshaped_NL_pred,block_length,width(reshaped_NL_pred)*length(reshaped_NL_pred)/block_length);

% throw our NaNs globally. so we're throwing out epochs where the data is incomplete
rm_this = isnan(sum(reshaped_LFP));
reshaped_LFP(:,rm_this) = [];
reshaped_PID(:,rm_this) = [];
reshaped_NL_pred(:,rm_this) = [];


% use the parameteric filter to project the stimulus, skipping the input nonlinearity 
[~, ~, ~, ~, K] = aNL4(reshaped_PID(:,1),p);
filtertime = 1e-3*(1:length(K));
NL_fp = NaN*reshaped_NL_pred;
time = 1e-3*(1:length(reshaped_LFP));
for i = 1:width(reshaped_LFP)
	NL_fp(:,i) = convolve(time,reshaped_PID(:,i),K,filtertime);	

	% remove mean from every trial
	NL_fp(:,i) = NL_fp(:,i) - nanmean(NL_fp(1e3:end,i));
end


subplot(4,4,13); hold on
[~,data_hi] = plotPieceWiseLinear(vectorise(NL_fp(1e3:5e3,:)),-vectorise(reshaped_LFP(1e3:5e3,:)),'Color',[1 0 0],'nbins',30);
[~,data_lo] = plotPieceWiseLinear(vectorise(NL_fp(6e3:9e3,:)),-vectorise(reshaped_LFP(6e3:9e3,:)),'Color',[0 0 1],'nbins',30);

% compute the model gain as a function of the standard deviation of the stimulus 

% correct the projected by a 1/S scaling
NL_fp_corrected = NL_fp;
for i = 1:width(NL_fp_corrected)
	NL_fp_corrected(1:5e3,i) = NLN_fp(1:5e3,i)/mean(reshaped_PID(1e3:4e3,i));
	NL_fp_corrected(5e3+1:end,i) = NLN_fp(5e3+1:end,i)/mean(reshaped_PID(6e3:9e3,i));
end

NL_fp_corrected = NL_fp_corrected*mean(reshaped_PID(:)); % overall correction to get the units right


% compute gains per trial on the corrected data
lo_gain = NaN(width(reshaped_PID),1);
hi_gain = NaN(width(reshaped_PID),1);
for i = 1:width(reshaped_PID)
	y = reshaped_NL_pred(1e3:4e3,i);
	x = NL_fp(1e3:4e3,i);
	try
		ff = fit(x(:),y(:),'poly1');
		hi_gain(i) = ff.p1;
	catch
	end

	y = reshaped_NL_pred(6e3:9e3,i);
	x = NL_fp(6e3:9e3,i);
	try
		ff = fit(x(:),y(:),'poly1');
		lo_gain(i) = ff.p1;
	catch
	end
end

lo_gain(lo_gain==0) = NaN;
hi_gain(hi_gain==0) = NaN;


% show gain as function of contrast
subplot(4,4,14); hold on
x = std(reshaped_PID(1e3:4e3,:));
plot(x,hi_gain,'r+')
x = std(reshaped_PID(6e3:9e3,:));
plot(x,lo_gain,'b+')
xlabel('\sigma_{Stimulus} (V)')
ylabel('ab3A ORN gain (Hz/V)')

%% Version Info
%
pFooter;





