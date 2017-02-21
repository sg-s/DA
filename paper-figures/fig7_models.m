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
p.     n = 1;


figure('outerposition',[0 0 800 999],'PaperUnits','points','PaperSize',[800 999]); hold on
clear ax
ax(1) = subplot(5,3,1); hold on
ax(2) = subplot(5,3,2); hold on
ax(3) = subplot(5,3,3); hold on



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
ylabel(ax(2),'Filter K (norm)')
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
ax(7) = subplot(5,3,7); hold on; cla
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
xlabel(ax(7),'Projected stimulus (V)')
ylabel(ax(7),'Model response (Hz)')
set(ax(7),'XLim',[0 2],'YLim',[0 60])

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

ax(8) = subplot(5,3,8); hold on; cla
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
	plot(ax(8),x,y,'+','Color',c(i,:));
end

% fit a power law with exponent -1
mean_stim = mean_stim(:);
g = MSGdata.NLN_gain(:);
options = fitoptions(fittype('power1'));
options.Lower = [-Inf -1];
options.Upper = [Inf -1];
cf = fit(mean_stim(~isnan(g)),g(~isnan(g)),'power1',options);
plot(ax(8),sort(mean_stim),cf(sort(mean_stim)),'r');
set(ax(8),'XScale','log','YScale','log','YLim',[10 300],'XLim',[.1 2.5])
xlabel(ax(8),'\mu_{Stimulus} (V)')
ylabel(ax(8),'Model gain (Hz/V)')

% compare model gains to actual gains
ax(9) = subplot(5,3,9); hold on; cla
for i = 1:max(MSGdata.paradigm)
	plot(ax(9),MSGdata.NLN_gain(MSGdata.paradigm == i),MSGdata.fA_gain(MSGdata.paradigm == i),'+','Color',c(i,:))
end
clear l
l = plot(ax(9),NaN,NaN,'k+');
legend(l,['r^2 = ' oval(rsquare(MSGdata.NLN_gain,MSGdata.fA_gain))],'Location','southeast')
xlabel(ax(9),'Model gain (Hz/V)')
ylabel(ax(9),'Observed gain (Hz/V)')
plot(ax(9),[0 250],[0 250],'k--')
set(ax(9),'XLim',[0 250],'YLim',[0 250])

% compute lags for each 
c = lines(5);
firing_color = c(4,:);


% compute the lags
lags = NaN*MSGdata.paradigm;
for i = 1:length(MSGdata.paradigm)
	S = MSGdata.PID(35e3:55e3,i); S = S - mean(S); S = S/std(S);
	R = MSGdata.NLN_pred(35e3:55e3,i); R = R - mean(R); R = R/std(R);
	lags(i) = finddelay(S,R);
end


y = NaN*(1:10);
ye = y;x = y;
for i = 1:10
	y(i) = nanmean(lags(MSGdata.paradigm == i));
	x(i) = nanmean(mean_stim(MSGdata.paradigm == i));
	ye(i) = nanstd(lags(MSGdata.paradigm == i));
end
subplot(5,3,15); hold on
set(gca,'YLim',[0 200])
errorbar(x,y,ye,'Color',firing_color)



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
ax(4) = subplot(5,3,4); hold on
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
legend(l,{'ab2A',['model r^2 = ' oval(rsquare(data(2).P(:,3),data(2).R(:,3)))]},'Location','northwest')

% show the output nonlinearity for ab3A
x = linspace(-.5,1,100);
y = x*p.C;
y(y<0) = 0;
clear l
l(1) = plot(ax(3),x,y,'r--');
l(2) = plot(ax(3),x,y*ff.p1,'r');
legend(l,{'ab3A','ab2A'},'Location','northwest')
set(ax(3),'XLim',[-.5 1 ],'YLim',[0 300])
xlabel(ax(3),'K \otimes a')
ylabel(ax(3),'Firing rate (Hz)')

xlabel(ax(4),'Time (s)')
ylabel(ax(4),'Firing rate (Hz)')

% show responses during the whiffs of approx same magnitude to show variation in response. 

ax(5) = subplot(5,3,5); hold on

show_these = [           2       27764
           2       28790
           2       59776
           3       58049];

for i = 2
	for j = 1:length(show_these)
		this_stim = show_these(j,1);
		this_loc = show_these(j,2);

		R = data(i).R(:,this_stim);

		a = this_loc - 300;
		z = this_loc+300;
		t = 1:length(R(a:z)); t = t - 300;
		plot(ax(5),t,R(a:z))

	end
end
ylabel(ax(5),'Firing rate (Hz)')
xlabel(ax(5),'Time since whiff (ms)')
set(ax(4),'XLim',[20 30])


% show r^2 for every whiff
ax(6) = subplot(5,3,6); hold on
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
l = plot(ax(6),x,y,'.','MarkerSize',20,'Color',[.5 .5 .5]);
legend(l,['r^2 = ' oval(rsquare(x,y))],'Location','southeast')
plot(ax(6),[0 300],[0 300],'k--')
xlabel(ax(6),' Model response (Hz)')
ylabel(ax(6),'ab2A response (Hz)')
title(ax(6),'Whiff-specific responses')
set(ax(6),'XLim',[0 300],'YLim',[0 300])

;;       ;;;;;;;; ;;;;;;;;     ;;     ;;  ;;;;;;   ;;;;;;   
;;       ;;       ;;     ;;    ;;;   ;;; ;;    ;; ;;    ;;  
;;       ;;       ;;     ;;    ;;;; ;;;; ;;       ;;        
;;       ;;;;;;   ;;;;;;;;     ;; ;;; ;;  ;;;;;;  ;;   ;;;; 
;;       ;;       ;;           ;;     ;;       ;; ;;    ;;  
;;       ;;       ;;           ;;     ;; ;;    ;; ;;    ;;  
;;;;;;;; ;;       ;;           ;;     ;;  ;;;;;;   ;;;;;;   


ax(15) = subplot(5,3,15); hold on

a = 35e3; z = 55e3; 
clear fd
for i = 1:length(MSGdata.paradigm)
	fd(i).stimulus = MSGdata.PID(a:z,i);
	R = MSGdata.raw_LFP(:,i);
	x = a:z;
	ff = fit(x(:),R(a:z),'poly1');
	R = R - ff(1:length(R));
	fd(i).response = R(a:z);
	fd(i).response(1:5e3) = NaN;
end

% generate responses using model
clear p
p.k0 = 0.1;
p.w = 0.6343;
p.B = 1.5984;
p.K_tau = 4.07;
p.output_offset = -0.508;
p.output_scale = -23.17;
p.n = 1;

% generate responses using the model 
if exist('.cache/asNL5_responses_to_MSG.mat','file') == 0
	
	LFP_pred = NaN*MSGdata.LFP;
	for i = 1:length(MSGdata.paradigm)
		textbar(i,length(MSGdata.paradigm))
		S = MSGdata.PID(:,i); S = S -min(S);
		LFP_pred(:,i) = asNL5(S,p);
	end
	save('.cache/asNL5_responses_to_MSG.mat','LFP_pred')
else
	load('.cache/asNL5_responses_to_MSG.mat')
end
MSGdata.LFP_pred = LFP_pred;
clear LFP_pred


% compute the lags
lags = NaN*MSGdata.paradigm;
for i = 1:length(MSGdata.paradigm)
	S = MSGdata.PID(35e3:55e3,i); S = S - mean(S); S = S/std(S);
	R = MSGdata.LFP_pred(35e3:55e3,i); R = R - mean(R); R = R/std(R);
	lags(i) = finddelay(S,R);
end

y = NaN*(1:10);
ye = y;x = y;
for i = 1:10
	y(i) = nanmean(lags(MSGdata.paradigm == i));
	x(i) = nanmean(mean_stim(MSGdata.paradigm == i));
	ye(i) = nanstd(lags(MSGdata.paradigm == i));
end
subplot(5,3,15); hold on
c = lines(5);
LFP_color = c(5,:);
errorbar(x,y,ye,'Color',LFP_color)
xlabel('\mu_{Stimulus} (V)')
ylabel('Lag (ms)')

;;       ;;;;;;;; ;;;;;;;;     ;;    ;;    ;;;    ;;;;;;;;     ;;;;;;  ;;;;;;;; ;;;; ;;     ;; 
;;       ;;       ;;     ;;    ;;;   ;;   ;; ;;      ;;       ;;    ;;    ;;     ;;  ;;;   ;;; 
;;       ;;       ;;     ;;    ;;;;  ;;  ;;   ;;     ;;       ;;          ;;     ;;  ;;;; ;;;; 
;;       ;;;;;;   ;;;;;;;;     ;; ;; ;; ;;     ;;    ;;        ;;;;;;     ;;     ;;  ;; ;;; ;; 
;;       ;;       ;;           ;;  ;;;; ;;;;;;;;;    ;;             ;;    ;;     ;;  ;;     ;; 
;;       ;;       ;;           ;;   ;;; ;;     ;;    ;;       ;;    ;;    ;;     ;;  ;;     ;; 
;;;;;;;; ;;       ;;           ;;    ;; ;;     ;;    ;;        ;;;;;;     ;;    ;;;; ;;     ;; 


% get all data 
cdata = consolidateData2(getPath(dataManager,'4608c42b12191b383c84fea52392ea97'));
[cdata, data] =  assembleScaledNatStim(cdata);


this_orn = 2;
clear fd
for i = 1:size(data(this_orn).X,2)
	S = data(this_orn).S(:,i); S = S - min(S);
	R = data(this_orn).X(:,i);
	fd(i).stimulus = S;
	fd(i).response = -R;
	fd(i).response(1:5e3) = NaN;
end

clear p
p.           k0 = 0.1000;
p.            w = 3.6250;
p.            B = 30.7500;
p.        K_tau = 21.7500;
p.output_offset = -0.0192;
p. output_scale = 24.7;
p.            n = 1;

% generate responses using this model 
for j = 1:size(data(2).X,2)
	data(2).XP(:,j) = -asNL5(data(2).S(:,j) - min(data(2).S(:,j)) ,p);
end




% show r^2 for every whiff
ax(13) = subplot(5,3,13); hold on
x = []; y = [];
for i = 2
	for j = 1:3
		S = data(i).S(:,j);
		P = data(i).XP(:,j);
		R = data(i).X(:,j);
		ws = whiffStatistics(S,R,R,300,'MinPeakProminence',max(S/1e2),'debug',false);
		y = [y; ws.peak_firing_rate];
		ws = whiffStatistics(S,P,P,300,'MinPeakProminence',max(S/1e2),'debug',false);
		x = [x; ws.peak_firing_rate];
	end
end
clear l

l = plot(ax(13),x,y,'.','MarkerSize',20,'Color',[.5 .5 .5]);
legend(l,['r^2 = ' oval(rsquare(x,y))],'Location','northwest')
xlabel(ax(13),'LFP model response (mV)')
ylabel(ax(13),'ab2 LFP response (mV)')
title(ax(13),'Whiff-specific responses')
plot(ax(13),[-18 0],[-18 0],'k--')
set(ax(13),'XLim',[-18 0],'YLim',[-18 0],'YDir','reverse','XDir','reverse')


% show context-dep. gain control
ax(14) = subplot(5,3,14); hold on
show_these = [2       27764
           2       28790
           2       59776
           3       58049];

for i = 2
	c = lines(length(show_these));
	for j = 1:length(show_these)
		this_stim = show_these(j,1);
		this_loc = show_these(j,2);

		P = data(i).XP(:,this_stim);

		a = this_loc - 300;
		z = this_loc+300;

		t = 1:length(P(a:z)); t = t - 300;
		plot(ax(14),t,P(a:z),'Color',c(j,:))

	end
end
ylabel(ax(14),'LFP (mV)')
xlabel(ax(14),'Time since whiff (ms)')
set(ax(14),'YDir','reverse')

% add the cartoon
ax_fig = subplot(5,3,10:12); hold on
axes(ax_fig)
o = imread('../images/lfp-model-cartoon.png');
imagesc(o);
axis ij
axis image
axis off

prettyFig('fs',12);

if being_published
	snapnow
	delete(gcf)
end



%% Version Info
%
pFooter;





