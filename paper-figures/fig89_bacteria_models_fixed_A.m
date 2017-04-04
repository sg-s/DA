%%
%  In this document, I fit a single adapting NLN model (aNLN4) to all the data.

pHeader;

% get data
clear MSGdata
MSGdata = consolidateData2(getPath(dataManager,'93ba5d68174e3df9f462a1fc48c581da'));
MSGdata = cleanMSGdata(MSGdata);

a = 35e3; z = 55e3; 

% make sure stimulus is always positive
for i = 1:size(MSGdata.PID,2)
	MSGdata.PID(:,i) = MSGdata.PID(:,i) - min(MSGdata.PID(:,i));
end
mean_stim = nanmean(MSGdata.PID(a:z,:));

fig1 = figure('outerposition',[0 0 1100 901],'PaperUnits','points','PaperSize',[1100 901]); hold on

;;       ;;;;;;;; ;;;;;;;;     ;;     ;;  ;;;;;;   ;;;;;;   
;;       ;;       ;;     ;;    ;;;   ;;; ;;    ;; ;;    ;;  
;;       ;;       ;;     ;;    ;;;; ;;;; ;;       ;;        
;;       ;;;;;;   ;;;;;;;;     ;; ;;; ;;  ;;;;;;  ;;   ;;;; 
;;       ;;       ;;           ;;     ;;       ;; ;;    ;;  
;;       ;;       ;;           ;;     ;; ;;    ;; ;;    ;;  
;;;;;;;; ;;       ;;           ;;     ;;  ;;;;;;   ;;;;;;   

figure(fig1)

% remove trends from the LFP 
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


% if you want to fit a model where A is fixed, use this:
clear p 
p.B = 2.0058;
p.e_L = 0;
p.K_1 = 0.09964;
p.K_2 = 10000;
p.K_tau = 2.0938;
p.n = 1;
p.output_scale = -22.066;
p.output_offset = 10.923;
p.Delta = 0;
p.Gamma = 0;
p.A = 10.68;

% generate responses using the model 
if exist('.cache/bacteriaModelX_fixed_A_responses_to_MSG.mat','file') == 0
	
	LFP_pred = NaN*MSGdata.LFP;
	for i = 1:length(MSGdata.paradigm)
		textbar(i,length(MSGdata.paradigm))
		S = MSGdata.PID(:,i); S = S -min(S);
		LFP_pred(:,i) = bacteriaModelX_fixedA(S,p);
	end
	save('.cache/bacteriaModelX_fixed_A_responses_to_MSG.mat','LFP_pred')
else
	load('.cache/bacteriaModelX_fixed_A_responses_to_MSG.mat')
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
ye = y; x = y;
for i = 1:10
	y(i) = nanmean(lags(MSGdata.paradigm == i));
	x(i) = nanmean(mean_stim(MSGdata.paradigm == i));
	ye(i) = nanstd(lags(MSGdata.paradigm == i));
end

subplot(3,3,6); hold on
c = lines(5);
LFP_color = c(5,:);
errorbar(x,y,ye,'Color',LFP_color)
xlabel('\mu_{Stimulus} (V)')
ylabel('Model lag (ms)')
box off


% gain vs. mean stimulus 
subplot(3,3,4);  hold on
c = parula(11);
gain = std(MSGdata.LFP_pred(a:z,:))./std(MSGdata.PID(a:z,:));
ff = fit(gain(:),MSGdata.LFP_gain(:),'poly1');
gain = ff(gain);
for i = 1:10
	plot(mean_stim(MSGdata.paradigm==i),gain(MSGdata.paradigm == i),'+','Color',c(i,:))
end
set(gca,'XScale','log','YScale','log')

text(.2, 3,'$\sim 1/s$','interpreter','latex','Color','r','FontSize',24)

% fit a power law with exponent -1
mean_stim = mean_stim(:);
g = gain(:);
options = fitoptions(fittype('power1'));
options.Lower = [-Inf -1];
options.Upper = [Inf -1];
cf = fit(mean_stim(~isnan(g)),g(~isnan(g)),'power1',options);
plot([min(mean_stim)/2 max(mean_stim)*2],cf([min(mean_stim)/2 max(mean_stim)*2]),'r');
set(gca,'XScale','log','YScale','log','YLim',[1 25],'XLim',[.1 2.5],'XTick',[.1 1],'YTick',[1 10])
xlabel('\mu_{Stimulus} (V)')
ylabel('Model gain (mV/V)')

% compare model gains to actual gains
subplot(3,3,5); hold on
for i = 1:max(MSGdata.paradigm)
	plot(gain(MSGdata.paradigm == i),MSGdata.LFP_gain(MSGdata.paradigm == i),'+','Color',c(i,:))
end
clear l
l = plot(NaN,NaN,'k+');
legend(l,['r^2 = ' oval(rsquare(gain,MSGdata.LFP_gain))],'Location','southeast')
xlabel('Model gain (mV/V)')
ylabel('LFP gain (mV/V)')
plot([0 250],[0 250],'k--')
set(gca,'XLim',[0 25],'YLim',[0 25])


;;       ;;;;;;;; ;;;;;;;;     ;;    ;;    ;;;    ;;;;;;;;     ;;;;;;  ;;;;;;;; ;;;; ;;     ;; 
;;       ;;       ;;     ;;    ;;;   ;;   ;; ;;      ;;       ;;    ;;    ;;     ;;  ;;;   ;;; 
;;       ;;       ;;     ;;    ;;;;  ;;  ;;   ;;     ;;       ;;          ;;     ;;  ;;;; ;;;; 
;;       ;;;;;;   ;;;;;;;;     ;; ;; ;; ;;     ;;    ;;        ;;;;;;     ;;     ;;  ;; ;;; ;; 
;;       ;;       ;;           ;;  ;;;; ;;;;;;;;;    ;;             ;;    ;;     ;;  ;;     ;; 
;;       ;;       ;;           ;;   ;;; ;;     ;;    ;;       ;;    ;;    ;;     ;;  ;;     ;; 
;;;;;;;; ;;       ;;           ;;    ;; ;;     ;;    ;;        ;;;;;;     ;;    ;;;; ;;     ;; 


% get all data 
cdata = consolidateData2(getPath(dataManager,'11f83bf78ccad7d44fee8e92dd65602f'));

% first remove min stimulus from every trace, and fix the LFP
for i = 1:size(cdata.PID,2)
	cdata.PID(:,i) = cdata.PID(:,i) - min(cdata.PID(:,i));
	cdata.LFP(:,i) = cdata.LFP(:,i) - mean(cdata.LFP(1:5e3,i));
end

cdata.LFP = cdata.LFP*10;

example_orn = 1;

S = mean(cdata.PID(:,cdata.orn == example_orn),2);
X = mean(cdata.LFP(:,cdata.orn == example_orn),2);
R = mean(cdata.fA(:,cdata.orn == example_orn),2);
time = 1e-3*(1:length(R));

S = S - min(S);


% used for fitting
clear fd
show_these = [28413 26669 64360 6112];
fd.stimulus = S;
fd.response = NaN*X;
for i = 1:length(show_these)
	a = show_these(i) - 200;
	z = show_these(i) + 200;
	fd.response(a:z) = X(a:z);
end

fd.response = X;
fd.response(1:5e3) = NaN;

clear p
p.B = 1.4126;
p.e_L = 0.9483;
p.w0 = 3000;
p.K_1 = 0.01;
p.K_2 = 1.5;
p.K_tau = 28.657;
p.n = 1;
p.output_offset = 8.889;
p.output_scale = -25.781;

% generate responses using this model 
XP = bacteriaModelX(S, p);

% show r^2 for every whiff
x = []; y = [];
ws = whiffStatistics(S,X,X,300,'MinPeakProminence',max(S/1e2),'debug',false);
y = [y; ws.peak_firing_rate];
ws = whiffStatistics(S,XP,XP,300,'MinPeakProminence',max(S/1e2),'debug',false);
x = [x; ws.peak_firing_rate];

clear l
subplot(3,3,9); hold on
l = plot(x,y,'.','MarkerSize',20,'Color',[.5 .5 .5]);
legend(l,['r^2 = ' oval(rsquare(x,y))],'Location','northwest')
xlabel('Model response (mV)')
ylabel('LFP response (mV)')
plot([-18 0],[-18 0],'k--')
set(gca,'XLim',[-18 0],'YLim',[-18 0],'YDir','reverse','XDir','reverse')


% show context-dep. gain control
subplot(3,3,8); hold on

show_these = [28413 26669 64360 6112];

for j = 1:length(show_these)
	this_loc = show_these(j);
	a = this_loc - 300;
	z = this_loc + 300;
	t = (1:length(XP(a:z))) - 300;

	plot(t,XP(a:z))

end
ylabel('Model response (mV)')
xlabel('Time since whiff (ms)')
set(gca,'YDir','reverse')

% also show the raw traces
subplot(3,3,7); hold on
plot(time,XP,'k')
plot(time,X,'r')
xlabel('Time (s)')
ylabel('\DeltaLFP (mV)')
set(gca,'XLim',[20 30],'YLim',[-25 5],'YDir','reverse')

clear th
th(1) = text(21, -20,'LFP','Color','k','FontSize',15);
th(2) = text(21, -17,'Model','Color','r','FontSize',15);

% add the cartoon
ax = subplot(3,3,1:3); hold on
ax.Position = [0 .68 1 .3];
axes(gca)
o = imread('../images/bacteria-model-eq.png');
imagesc(o);
axis ij
axis image
axis off

prettyFig('fs',15);

[~,lh] = labelFigure('x_offset',0);
lh(1).Position = [.07 .93 .02 .03];

if being_published
	snapnow
	delete(gcf)
end

return

;;;;;;;; ;;;; ;;;;;;;;  ;;;; ;;    ;;  ;;;;;;      ;;;;;;;;     ;;;    ;;;;;;;; ;;;;;;;; 
;;        ;;  ;;     ;;  ;;  ;;;   ;; ;;    ;;     ;;     ;;   ;; ;;      ;;    ;;       
;;        ;;  ;;     ;;  ;;  ;;;;  ;; ;;           ;;     ;;  ;;   ;;     ;;    ;;       
;;;;;;    ;;  ;;;;;;;;   ;;  ;; ;; ;; ;;   ;;;;    ;;;;;;;;  ;;     ;;    ;;    ;;;;;;   
;;        ;;  ;;   ;;    ;;  ;;  ;;;; ;;    ;;     ;;   ;;   ;;;;;;;;;    ;;    ;;       
;;        ;;  ;;    ;;   ;;  ;;   ;;; ;;    ;;     ;;    ;;  ;;     ;;    ;;    ;;       
;;       ;;;; ;;     ;; ;;;; ;;    ;;  ;;;;;;      ;;     ;; ;;     ;;    ;;    ;;;;;;;; 


fig2 = figure('outerposition',[0 0 1100 901],'PaperUnits','points','PaperSize',[1100 901]); hold on
clear ax

clear p
p.     B =  0.8289;
p.tau_y1 =  27;
p.tau_y2 =  96;
p.     A =  0.7000;
p.     n =  1;
p.     C =  178;
p.   K_1 =  1.0000e-06;
p.   e_L =  11.5129;
p.   K_2 =  1000000;

% first, show 3 nonlinearities from the Gaussians with the right colours 
ax(1) = subplot(3,3,1); hold on
c = parula(11);
for i = [2 6 10]
	this_paradigm = find(MSGdata.paradigm == i,1,'first');
	S = MSGdata.PID(:,this_paradigm);
	[~, a, ~, e0] = bacteriaModelF(S,p);

	x = logspace(-2,2,1e3);
	Shat = (1 + x/p.K_2)./(1 + x/p.K_1);
	E = exp(mean(e0(35e3:55e3)) + log(Shat));
	a_mean = 1./(1 + E);
	plot(x,a_mean,'Color',c(i,:))

	plot(S(35e3:55e3),a(35e3:55e3),'.','Color',c(i,:))
end
set(gca,'XScale','log','XLim',[1e-2 10],'YAxisLocation','right')
xlabel('Stimulus (V)')
ylabel('a')
axis square

% show the filter 
q.n = p.n;
q.tau1 = p.tau_y1;
q.tau2 = p.tau_y2;
q.A = p.A;
K = filter_gamma2(1:1e3,q);
ax(2) = subplot(3,3,2); hold on
plot(K/norm(K),'k')
set(gca,'XLim',[0 400],'YLim',[-.05 .15])
xlabel('Model lag (ms)')
ylabel('Filter K (norm)')
axis square

% show the output nonlinearity 
x = linspace(-1,2,1e3);
y = x*p.C; y(y<0) = 0;
ax(3) = subplot(3,3,3); hold on
plot(x,y,'k')
xlabel('K \otimes a')
ylabel('Model response (Hz)')
set(gca,'XLim',[-1 2],'YLim',[-5 250])
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
if exist('.cache/bacteriaModelF_responses_to_MSG.mat','file') == 0
	NLN_pred = NaN*MSGdata.PID;
	for i = 1:length(MSGdata.paradigm)
		S = MSGdata.PID(:,i) - min(MSGdata.PID(:,i));
		NLN_pred(:,i) = bacteriaModelF(S,p);
	end
	save('.cache/bacteriaModelF_responses_to_MSG.mat','NLN_pred')
else
	load('.cache/bacteriaModelF_responses_to_MSG.mat')
end
MSGdata.NLN_pred = NLN_pred;
clear NLN_pred

% back out new linear filters for this
time = 1e-3*(1:length(MSGdata.PID));
if exist('.cache/bacteriaModelF_MSG_linear_prediction.mat','file') == 0
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
	save('.cache/bacteriaModelF_MSG_linear_prediction.mat','NLN_fp')
else
	load('.cache/bacteriaModelF_MSG_linear_prediction.mat')
end
MSGdata.NLN_fp = NLN_fp;
clear NLN_fp

% recreate the I/O plot 
ax(4) = subplot(3,3,4); hold on
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
ylabel('Model response (Hz)')
set(gca,'XLim',[0 2],'YLim',[0 60])


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

% correct for some trivial scaling 
x = MSGdata.NLN_gain(:);
y = MSGdata.fA_gain(:);
rm_this = isnan(x) | isnan(y);
ff = fit(x(~rm_this),y(~rm_this),'poly1');
MSGdata.NLN_gain = ff(MSGdata.NLN_gain);

% show gain changes -- gain vs. mean stimulus
ax(5) = subplot(3,3,5); hold on
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
ylabel('Model gain (Hz/V)')

% compare model gains to actual gains
ax(6) = subplot(3,3,6); hold on
for i = 1:max(MSGdata.paradigm)
	plot(MSGdata.NLN_gain(MSGdata.paradigm == i),MSGdata.fA_gain(MSGdata.paradigm == i),'+','Color',c(i,:))
end
clear l
l = plot(NaN,NaN,'k+');
legend(l,['r^2 = ' oval(rsquare(MSGdata.NLN_gain,MSGdata.fA_gain))],'Location','southeast')
xlabel('Model gain (Hz/V)')
ylabel('Observed gain (Hz/V)')
plot([0 250],[0 250],'k--')
set(gca,'XLim',[0 250],'YLim',[0 250])

% % compute lags for each 
% c = lines(5);
% firing_color = c(4,:);


% % compute the lags
% lags = NaN*MSGdata.paradigm;
% for i = 1:length(MSGdata.paradigm)
% 	S = MSGdata.PID(35e3:55e3,i); S = S - mean(S); S = S/std(S);
% 	R = MSGdata.NLN_pred(35e3:55e3,i); R = R - mean(R); R = R/std(R);
% 	lags(i) = finddelay(S,R);
% end


% y = NaN*(1:10);
% ye = y;x = y;
% for i = 1:10
% 	y(i) = nanmean(lags(MSGdata.paradigm == i));
% 	x(i) = nanmean(mean_stim(MSGdata.paradigm == i));
% 	ye(i) = nanstd(lags(MSGdata.paradigm == i));
% end
% subplot(5,3,15); hold on
% set(gca,'YLim',[0 200])
% errorbar(x,y,ye,'Color',firing_color)

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

clear p
p.     B =  0.8289;
p.tau_y1 =  27;
p.tau_y2 =  96;
p.     A =  0.7000;
p.     n =  1;
p.     C =  178;
p.   K_1 =  1.0000e-06;
p.   e_L =  11.5129;
p.   K_2 =  1000000;

% show the naturalistic stimulus fits
cdata = consolidateData2(getPath(dataManager,'11f83bf78ccad7d44fee8e92dd65602f'));

% first remove min stimulus from every trace, and fix the LFP
for i = 1:size(cdata.PID,2)
	cdata.PID(:,i) = cdata.PID(:,i) - min(cdata.PID(:,i));
end

example_orn = 1;

S = mean(cdata.PID(:,cdata.orn == example_orn),2);
R = mean(cdata.fA(:,cdata.orn == example_orn),2);
time = 1e-3*(1:length(R));
S = S - min(S);


% show the time series of the natualistic stimulus
ax(7) = subplot(3,3,7); hold on
clear l
l(1) = plot(ax(7),time,R,'k');

% generate responses and plot them too
P = bacteriaModelF(S - min(S) ,p);

% fix some trivial scaling
ff = fit(P(:),R,'poly1');
P = ff(P); P(P<0) = 0;

l(2) = plot(ax(7),time,P,'r');


xlabel(ax(7),'Time (s)')
ylabel(ax(7),'Firing rate (Hz)')
set(ax(7),'XLim',[20 30])

% show responses during the whiffs of approx same magnitude to show variation in response. 

ax(8) = subplot(3,3,8); hold on

show_these = [28413 26669 64360 6112];


for j = 1:length(show_these)
	this_loc = show_these(j);

	tR = mean(P,2);

	a = this_loc - 300;
	z = this_loc + 300;

	t = (1:length(tR(a:z))) - 300;

	plot(ax(8),t,tR(a:z))

end

ylabel(ax(8),'Firing rate (Hz)')
xlabel(ax(8),'Time since whiff (ms)')


% show r^2 for every whiff

x = []; y = [];

ws = whiffStatistics(S,R,R,300,'MinPeakProminence',max(S/1e2),'debug',false);
y = [y; ws.peak_firing_rate];
ws = whiffStatistics(S,P,P,300,'MinPeakProminence',max(S/1e2),'debug',false);
x = [x; ws.peak_firing_rate];

clear l
ax(9) = subplot(3,3,9); hold on
l = plot(ax(9),x,y,'.','MarkerSize',20,'Color',[.5 .5 .5]);
legend(l,['r^2 = ' oval(rsquare(x,y))],'Location','southeast')
plot(ax(9),[0 300],[0 300],'k--')
xlabel(ax(9),' Model response (Hz)')
ylabel(ax(9),'ORN response (Hz)')
set(ax(9),'XLim',[50 200],'YLim',[50 200])


prettyFig('fs',15);

% add labels
[~,lh] = labelFigure('x_offset',0,'y_offset',0.01);


if being_published
	snapnow
	delete(gcf)
end



%% Version Info
%
pFooter;





