% fig1_nat_stim_example
% now uses ab2A and ab3A data 

%%
% This script makes fig 1 for the eLife resubmission. The point of this figure is to show our naturalistic stimulus data, and the features of the ORN response we see in it. 

pHeader;


% make the figure and all subplots 
fig1 = figure('outerposition',[0 0 1400 777],'PaperUnits','points','PaperSize',[1400 777]); hold on
clear ax

% time series for ab3A
ax.ab3A_S = subplot(6,6,1:3); hold on
ax.ab3A_X = subplot(6,6,7:9); hold on
ax.ab3A_R = subplot(6,6,13:15); hold on

% time series for ab2A
ax.ab2A_S = subplot(6,6,19:21); hold on
ax.ab2A_X = subplot(6,6,25:27); hold on
ax.ab2A_R = subplot(6,6,31:33); hold on

% filters + projections for ab2A 
ax.X_proj = subplot(2,4,3); hold on
ax.R_proj = subplot(2,4,4); hold on

% dose response for ab2A
ax.ab2A_drX = subplot(2,4,7); hold on
ax.ab2A_drR = subplot(2,4,8); hold on

fig2 = figure('outerposition',[0 0 1400 777],'PaperUnits','points','PaperSize',[1400 777]); hold on

% zooms into time series, one for ab3A, one for ab2A 
ax.ab3A_SZ = subplot(6,6,1); hold on
ax.ab3A_XZ = subplot(6,6,7); hold on
ax.ab3A_RZ = subplot(6,6,13); hold on
ax.ab2A_SZ = subplot(6,6,19); hold on
ax.ab2A_XZ = subplot(6,6,25); hold on
ax.ab2A_RZ = subplot(6,6,31); hold on


% add insets showing that whiffs don't change, but responses do
ax.ab3A_SZi = subplot(6,6,2); hold on
ax.ab3A_SZi2 = subplot(6,6,3); hold on
ax.ab3A_XZi = subplot(6,6,9); hold on
ax.ab3A_RZi = subplot(6,6,15); hold on


ax.ab2A_SZi = subplot(6,6,20); hold on
ax.ab2A_SZi2 = subplot(6,6,21); hold on
ax.ab2A_XZi = subplot(6,6,27); hold on
ax.ab2A_RZi = subplot(6,6,33); hold on

ax.deviations_X_300ms = subplot(2,4,3); hold on
ax.deviations_R_300ms = subplot(2,4,4); hold on

ax.deviations_X_0ms = axes;
ax.deviations_X_0ms.Position = [.6 .8 .1 .1];

ax.deviations_R_0ms = axes;
ax.deviations_R_0ms.Position = [.8 .8 .1 .1];

ax.deviations_prev_whiff_X = subplot(2,4,7); hold on
ax.deviations_prev_whiff_R = subplot(2,4,8); hold on

   ;;;    ;;;;;;;;   ;;;;;;;     ;;;    
  ;; ;;   ;;     ;; ;;     ;;   ;; ;;   
 ;;   ;;  ;;     ;;        ;;  ;;   ;;  
;;     ;; ;;;;;;;;   ;;;;;;;  ;;     ;; 
;;;;;;;;; ;;     ;; ;;        ;;;;;;;;; 
;;     ;; ;;     ;; ;;        ;;     ;; 
;;     ;; ;;;;;;;;  ;;;;;;;;; ;;     ;; 


% get the data

cdata = consolidateData2(getPath(dataManager,'4608c42b12191b383c84fea52392ea97'));
[cdata, data] =  assembleScaledNatStim(cdata);
time = 1e-3*(1:length(data(1).S));

% plot the data
plot(ax.ab2A_S,time,data(2).S(:,1),'k')
plot(ax.ab2A_X,time,data(2).X(:,1),'k')
plot(ax.ab2A_R,time,data(2).R(:,1),'k')

% show whiff statistics 
i = 2;
clear all_x all_y
all_x = []; all_y = [];
for j = 1:size(data(i).S,2)
	S = data(i).S(:,j);
	X = data(i).X(:,j);
	R = data(i).R(:,j);
	ws = whiffStatistics(S,X,R,300,'MinPeakProminence',max(S/1e2),'debug',false);
	all_x =  [all_x(:); ws.stim_peaks(:)];
	all_y = [all_y(:); -ws.peak_LFP(:)];
	plot(ax.ab2A_drX,ws.stim_peaks,ws.peak_LFP,'.','MarkerSize',20,'Color','k')
	plot(ax.ab2A_drR,ws.stim_peaks,ws.peak_firing_rate,'.','MarkerSize',20,'Color','k')
end

% show context-dependent variation in whiff response
s_range = [.8 1.2];

show_these = [2    28790
           2       27764
           2       59776
           3       58049];

for i = 2
	for j = 1:length(show_these)
		this_stim = show_these(j,1);
		this_loc = show_these(j,2);

		S = data(i).S(:,this_stim);
		X = data(i).X(:,this_stim);
		R = data(i).R(:,this_stim);

		a = this_loc - 300;
		z = this_loc+300;

		t = (1:length(S(a:z))) - 300;

		plot(ax.ab2A_SZ,t,S(a:z))
		plot(ax.ab2A_XZ,t,X(a:z))
		plot(ax.ab2A_RZ,t,R(a:z))

	end
end

% show some statistics about these chosen whiffs
S_int = NaN*(1:4);
S_int_err = NaN*(1:4);
X_int = NaN*(1:4);
X_int_err = NaN*(1:4);
R_int = NaN*(1:4);
R_int_err = NaN*(1:4);

S_before = NaN*(1:4);
S_before_err = NaN*(1:4);

for i = 1:length(show_these)
	this_paradigm = show_these(i,1); 

	PID = data(2).S(:,this_paradigm);
	LFP = data(2).X(:,this_paradigm);
	fA =  data(2).R(:,this_paradigm);

	temp = max(PID(show_these(i,2)-300:show_these(i,2)-100,:));
	S_before(i) = nanmean(temp);
	S_before_err(i) = nanstd(temp);

	temp = nanmean(PID(show_these(i,2):show_these(i,2)+10,:));
	S_int(i) = nanmean(temp);
	S_int_err(i) = nanstd(temp);

	temp = nanmean(LFP(show_these(i,2)+90:show_these(i,2)+100,:));
	X_int(i) = nanmean(temp);
	X_int_err(i) = nanstd(temp);

	temp = nanmean(fA(show_these(i,2)+60:show_these(i,2)+80,:));
	R_int(i) = nanmean(temp);
	R_int_err(i) = nanstd(temp);

end

cla(ax.ab2A_SZi)
cla(ax.ab2A_XZi)
cla(ax.ab2A_RZi)
c = lines(10);

% show whiff statistics
for i = 1:length(show_these)
	superbar(ax.ab2A_SZi2,(i),S_int(i),'E',S_int_err(i),'BarFaceColor',c(i,:),'ErrorbarColor',c(i,:))
	superbar(ax.ab2A_XZi,(i),X_int(i),'E',X_int_err(i),'BarFaceColor',c(i,:),'ErrorbarColor',c(i,:))
	superbar(ax.ab2A_RZi,(i),R_int(i),'E',R_int_err(i),'BarFaceColor',c(i,:),'ErrorbarColor',c(i,:))
end

% also show the statistics of the whiff before
for i = 1:length(show_these)
	superbar(ax.ab2A_SZi,(i),(S_before(i)),'E',abs((S_before_err(i))),'BarFaceColor',c(i,:),'ErrorbarColor',c(i,:))

end

set(ax.ab2A_SZi,'YLim',[0 1.5])
set(ax.ab2A_SZi2,'YLim',[0 1.5])


% show response vs. projected stimulus
K = fitFilter2Data(data(2).S(:,1),data(2).X(:,1),'offset',200); 
K = K(100:end-100);
filtertime = (1:length(K)) - 100;
fp = convolve(1:length(data(2).S(:)),data(2).S(:),K,filtertime);
plot(ax.X_proj,fp,data(2).X(:),'k');
set(ax.X_proj,'YDir','reverse','XDir','reverse','XLim',[-4 0],'YLim',[-30 0])
xlabel(ax.X_proj,'Projected stimulus (V)')
ylabel(ax.X_proj,'ab2 LFP (mV)')
l = plot(ax.X_proj,NaN,NaN,'Color',[ 1 1 1]);
legend(l,['r^2 = ' oval(rsquare(fp,data(2).X(:)))],'Location','southeast');

K = fitFilter2Data(data(2).S(:,1),data(2).R(:,1),'offset',200); 
K = K(100:end-100);
filtertime = (1:length(K)) - 100;
fp = convolve(1:length(data(2).S(:)),data(2).S(:),K,filtertime);
plot(ax.R_proj,fp,data(2).R(:),'k');
set(ax.R_proj,'XLim',[0 5],'YLim',[0 300])
xlabel(ax.R_proj,'Projected stimulus (V)')
ylabel(ax.R_proj,'ab2A firing rate (Hz)')
l = plot(ax.R_proj,NaN,NaN,'Color',[ 1 1 1]);
legend(l,['r^2 = ' oval(rsquare(fp,data(2).R(:)))],'Location','southeast');

;;;;;;;;  ;;;;;;;; ;;     ;; ;;;;    ;;;    ;;;;;;;; ;;;;  ;;;;;;;  ;;    ;; 
;;     ;; ;;       ;;     ;;  ;;    ;; ;;      ;;     ;;  ;;     ;; ;;;   ;; 
;;     ;; ;;       ;;     ;;  ;;   ;;   ;;     ;;     ;;  ;;     ;; ;;;;  ;; 
;;     ;; ;;;;;;   ;;     ;;  ;;  ;;     ;;    ;;     ;;  ;;     ;; ;; ;; ;; 
;;     ;; ;;        ;;   ;;   ;;  ;;;;;;;;;    ;;     ;;  ;;     ;; ;;  ;;;; 
;;     ;; ;;         ;; ;;    ;;  ;;     ;;    ;;     ;;  ;;     ;; ;;   ;;; 
;;;;;;;;  ;;;;;;;;    ;;;    ;;;; ;;     ;;    ;;    ;;;;  ;;;;;;;  ;;    ;; 

   ;;;    ;;    ;;    ;;;    ;;       ;;    ;;  ;;;;;;  ;;;;  ;;;;;;  
  ;; ;;   ;;;   ;;   ;; ;;   ;;        ;;  ;;  ;;    ;;  ;;  ;;    ;; 
 ;;   ;;  ;;;;  ;;  ;;   ;;  ;;         ;;;;   ;;        ;;  ;;       
;;     ;; ;; ;; ;; ;;     ;; ;;          ;;     ;;;;;;   ;;   ;;;;;;  
;;;;;;;;; ;;  ;;;; ;;;;;;;;; ;;          ;;          ;;  ;;        ;; 
;;     ;; ;;   ;;; ;;     ;; ;;          ;;    ;;    ;;  ;;  ;;    ;; 
;;     ;; ;;    ;; ;;     ;; ;;;;;;;;    ;;     ;;;;;;  ;;;;  ;;;;;;  


S = data(2).S; 
X = data(2).X;
R = data(2).R;
S = S(:);

clear ws
for j = 1:3
	ws(j) = whiffStatistics(data(2).S(:,j),data(2).X(:,j),data(2).R(:,j),300,'MinPeakProminence',max(data(2).S(:,j)/100));
end

for i = 1:2
	ws(i+1).stim_peak_loc = ws(i+1).stim_peak_loc + 70e3*i;
end

stim_peaks = vertcat(ws.stim_peaks);
stim_peak_loc = vertcat(ws.stim_peak_loc);
R_peak = vertcat(ws.peak_firing_rate);
X_peak = vertcat(ws.peak_LFP);

% compute the "median" for each whiff, and the deviations from the median for each whiff 
M_X = NaN*stim_peaks;
M_R = NaN*stim_peaks;
bin_width = .2;
for i = 1:length(stim_peaks)
	other_X = X_peak(stim_peaks > stim_peaks(i)*(1-bin_width) & stim_peaks < stim_peaks(i)*(1+bin_width));
	other_R = R_peak(stim_peaks > stim_peaks(i)*(1-bin_width) & stim_peaks < stim_peaks(i)*(1+bin_width));
	if length(other_R) > 3
		M_R(i) = median(other_R);
		M_X(i) = -median(other_X);
	end
end

D_X = (-X_peak - M_X)./M_X;
D_X(D_X==0) = NaN;
D_R = (R_peak - M_R)./M_R;
D_R(D_R==0) = NaN;


%  plot deviations as a function of whiff intensity to show no correlation 
l = plot(ax.deviations_X_0ms,stim_peaks,D_X,'.','Color',[.5 .5 .5],'MarkerSize',24);
[rho,p]=spear(stim_peaks,D_X);
lh.x_no_corr = legend(l,['\rho = ' oval(rho), ', p = ' oval(p)]);
set(ax.deviations_X_0ms,'XScale','log')
xlabel(ax.deviations_X_0ms,'Whiff amplitude (V)')
ylabel(ax.deviations_X_0ms,'D_{LFP}')

l = plot(ax.deviations_R_0ms,stim_peaks,D_R,'.','Color',[.5 .5 .5],'MarkerSize',24);
[rho,p]=spear(stim_peaks,D_R);
lh.r_no_corr = legend(l,['\rho = ' oval(rho), ', p = ' oval(p)]);
set(ax.deviations_R_0ms,'XScale','log')
xlabel(ax.deviations_R_0ms,'Whiff amplitude (V)')
ylabel(ax.deviations_R_0ms,'D_{F}')


% plot deviations as a function of the preceding stimulus in 300ms
Shat = computeSmoothedStimulus(S,300)';
l = plot(ax.deviations_X_300ms,Shat(stim_peak_loc),D_X,'.','Color','k','MarkerSize',20);
[rho,p]=spear(Shat(stim_peak_loc),D_X);
lh.x_300 = legend(l,['\rho = ' oval(rho), ', p = ' oval(p)]);
xlabel(ax.deviations_X_300ms,'\mu_{stimulus} in preceding 300ms (V)')
ylabel(ax.deviations_X_300ms,'Fractional LFP deviations D_{LFP}')

Shat = computeSmoothedStimulus(S,300)';
l = plot(ax.deviations_R_300ms,Shat(stim_peak_loc),D_R,'.','Color','k','MarkerSize',20);
[rho,p]=spear(Shat(stim_peak_loc),D_R);
lh.r_300 = legend(l,['\rho = ' oval(rho), ', p < 10^{-5}']);
xlabel(ax.deviations_R_300ms,'\mu_{stimulus} in preceding 300ms (V)')
ylabel(ax.deviations_R_300ms,'Fractional firing rate deviations D_{F}')

% plot deviations as a function of previous whiff 
T_before = stim_peak_loc - circshift(stim_peak_loc,1); % the time to the preceding whiff
S_before = circshift(stim_peaks,1); % the peak of the preceding whiff
X_low = [(S_before(D_X<0)) T_before(D_X<0)];
X_high = [(S_before(D_X>0)) T_before(D_X>0)];
R_low = [(S_before(D_R<0)) T_before(D_R<0)];
R_high = [(S_before(D_R>0)) T_before(D_R>0)];

clear l
l(1) = plot(ax.deviations_prev_whiff_X,X_high(:,1),X_high(:,2),'r.','MarkerSize',48);
l(2) = plot(ax.deviations_prev_whiff_X,X_low(:,1),X_low(:,2),'b.','MarkerSize',48);
[~, p] = kstest_2s_2d(X_low, X_high);
legend(l,{'D_{LFP} > 0','D_{LFP} < 0'})
set(ax.deviations_prev_whiff_X,'XScale','log','YScale','log')
xlabel(ax.deviations_prev_whiff_X,'Amplitude of previous whiff (V)')
ylabel(ax.deviations_prev_whiff_X,'Time since previous whiff (ms)')

% display some statistics
disp('2D KS test for LFP: ')
disp(['p = ' oval(p)])
disp('KS test using whiff amplitudes only:')
[~,p] = kstest2(X_high(:,1),X_low(:,1));
disp(['p = ' oval(p)])
disp('KS test using whiff times only:')
[~,p] = kstest2(X_high(:,2),X_low(:,2));
disp(['p = ' oval(p)])


clear l
l(1) = plot(ax.deviations_prev_whiff_R,R_high(:,1),R_high(:,2),'r.','MarkerSize',48);
l(2) = plot(ax.deviations_prev_whiff_R,R_low(:,1),R_low(:,2),'b.','MarkerSize',48);
[~, p] = kstest_2s_2d(R_low, R_high);
legend(l,{'D_F > 0','D_F < 0'})
set(ax.deviations_prev_whiff_R,'XScale','log','YScale','log')
xlabel(ax.deviations_prev_whiff_R,'Amplitude of previous whiff (V)')
ylabel(ax.deviations_prev_whiff_R,'Time since previous whiff (ms)')

% display some statistics
disp('2D KS test for firing rate: ')
disp(['p = ' oval(p)])
disp('KS test using whiff amplitudes only:')
[~,p] = kstest2(R_high(:,1),R_low(:,1));
disp(['p = ' oval(p)])
disp('KS test using whiff times only:')
[~,p] = kstest2(R_high(:,2),R_low(:,2));
disp(['p = ' oval(p)])



   ;;;    ;;;;;;;;   ;;;;;;;     ;;;    
  ;; ;;   ;;     ;; ;;     ;;   ;; ;;   
 ;;   ;;  ;;     ;;        ;;  ;;   ;;  
;;     ;; ;;;;;;;;   ;;;;;;;  ;;     ;; 
;;;;;;;;; ;;     ;;        ;; ;;;;;;;;; 
;;     ;; ;;     ;; ;;     ;; ;;     ;; 
;;     ;; ;;;;;;;;   ;;;;;;;  ;;     ;; 


% now do the ab3A analysis
cdata = consolidateData2(getPath(dataManager,'11f83bf78ccad7d44fee8e92dd65602f'));


% first remove min stimulus from every trace, and fix the LFP
for i = 1:size(cdata.PID,2)
	cdata.PID(:,i) = cdata.PID(:,i) - min(cdata.PID(:,i));
	cdata.LFP(:,i) = cdata.LFP(:,i) - mean(cdata.LFP(1:5e3,i));
end

cdata.LFP = cdata.LFP*10;

example_orn = 1;

PID = cdata.PID(:,cdata.orn == example_orn);
LFP = cdata.LFP(:,cdata.orn == example_orn);
fA = cdata.fA(:,cdata.orn == example_orn);

S = mean(cdata.PID(:,cdata.orn == example_orn),2);
X = mean(cdata.LFP(:,cdata.orn == example_orn),2);
R = mean(cdata.fA(:,cdata.orn == example_orn),2);

rm_this = isnan(sum(S)) | max(R) == 0;
S(:,rm_this) = []; X(:,rm_this) = []; R(:,rm_this) = [];

% plot the data
plot(ax.ab3A_S,time,mean(S,2),'k')
plot(ax.ab3A_X,time,mean(X,2),'k')
plot(ax.ab3A_R,time,mean(R,2),'k')


clear all_x all_y
all_x = []; all_y = [];




% now show fast gain control
% show context-dependent variation in whiff response

% show_these = [21103 64353 64869];
show_these = [28413 26669 64360 6112];

cla(ax.ab3A_SZ)
cla(ax.ab3A_XZ)
cla(ax.ab3A_RZ)


for j = 1:length(show_these)
	this_loc = show_these(j);

	tS = mean(S,2);
	tX = mean(X,2);
	tR = mean(R,2);

	a = this_loc - 300;
	z = this_loc + 300;

	t = (1:length(tS(a:z))) - 300;

	plot(ax.ab3A_SZ,t,tS(a:z))
	plot(ax.ab3A_XZ,t,tX(a:z))
	plot(ax.ab3A_RZ,t,tR(a:z))

end


% show some statistics about these chosen whiffs
S_int = NaN*show_these;
S_int_err = NaN*show_these;
X_int = NaN*show_these;
X_int_err = NaN*show_these;
R_int = NaN*show_these;
R_int_err = NaN*show_these;

S_before = NaN*(1:4);
S_before_err = NaN*(1:4);

for i = 1:length(show_these)
	temp = nanmean(PID(show_these(i):show_these(i)+10,:));
	S_int(i) = nanmean(temp);
	S_int_err(i) = nanstd(temp)/sqrt(length(temp));

	temp = max(PID(show_these(i)-300:show_these(i)-100,:));
	S_before(i) = nanmean(temp);
	S_before_err(i) = nanstd(temp);


	temp = nanmean(LFP(show_these(i)+90:show_these(i)+100,:));
	X_int(i) = nanmean(temp);
	X_int_err(i) = nanstd(temp)/sqrt(length(temp));

	temp = nanmean(fA(show_these(i)+90:show_these(i)+100,:));
	R_int(i) = nanmean(temp);
	R_int_err(i) = nanstd(temp)/sqrt(length(temp));

end

cla(ax.ab3A_SZi)
cla(ax.ab3A_XZi)
cla(ax.ab3A_RZi)
c = lines(10);


for i = 1:length(show_these)
	superbar(ax.ab3A_SZi2,i,S_int(i),'E',S_int_err(i),'BarFaceColor',c(i,:),'ErrorbarColor',c(i,:))
	superbar(ax.ab3A_SZi,i,S_before(i),'E',S_before_err(i),'BarFaceColor',c(i,:),'ErrorbarColor',c(i,:))
	superbar(ax.ab3A_XZi,i,X_int(i),'E',X_int_err(i),'BarFaceColor',c(i,:),'ErrorbarColor',c(i,:))
	superbar(ax.ab3A_RZi,i,R_int(i),'E',R_int_err(i),'BarFaceColor',c(i,:),'ErrorbarColor',c(i,:))

end

 ;;;;;;   ;;;;;;;   ;;;;;;  ;;     ;; ;;;;;;;; ;;;;;;;; ;;;;  ;;;;;;   ;;;;;;  
;;    ;; ;;     ;; ;;    ;; ;;;   ;;; ;;          ;;     ;;  ;;    ;; ;;    ;; 
;;       ;;     ;; ;;       ;;;; ;;;; ;;          ;;     ;;  ;;       ;;       
;;       ;;     ;;  ;;;;;;  ;; ;;; ;; ;;;;;;      ;;     ;;  ;;        ;;;;;;  
;;       ;;     ;;       ;; ;;     ;; ;;          ;;     ;;  ;;             ;; 
;;    ;; ;;     ;; ;;    ;; ;;     ;; ;;          ;;     ;;  ;;    ;; ;;    ;; 
 ;;;;;;   ;;;;;;;   ;;;;;;  ;;     ;; ;;;;;;;;    ;;    ;;;;  ;;;;;;   ;;;;;;  

% clean up fig 1
axis(ax.X_proj,'square')
axis(ax.R_proj,'square')
axis(ax.ab2A_drX,'square')
axis(ax.ab2A_drR,'square')

set(ax.ab2A_drX,'XScale','log','XLim',[1e-2 1e1],'YDir','reverse','XTick',[1e-2 1e-1 1e0 1e1 1e2])
xlabel(ax.ab2A_drR,'Whiff amplitude (V)')
ylabel(ax.ab2A_drR,'ab2A firing rate (Hz)')
xlabel(ax.ab2A_drX,'Whiff amplitude (V)')
ylabel(ax.ab2A_drX,'ab2 LFP rate (mV)')
set(ax.ab2A_drR,'XScale','log','XLim',[1e-2 1e1],'XTick',[1e-2 1e-1 1e0 1e1 1e2])

ylabel(ax.ab3A_S,['Stimulus' char(10) '(V)'])
set(ax.ab3A_S,'XLim',[0 70],'XTick',[])
ylabel(ax.ab3A_X,['\DeltaLFP' char(10) '(mV)'])
set(ax.ab3A_X,'XLim',[0 70],'YLim',[-20 2])
set(ax.ab3A_X,'YDir','reverse','XTick',[])
ylabel(ax.ab3A_R,['Firing' char(10) 'rate (Hz)'])
set(ax.ab3A_R,'XLim',[0 70])
ylabel(ax.ab2A_S,['Stimulus' char(10) '(V)'])
set(ax.ab2A_S,'XLim',[0 70],'XTick',[])
ylabel(ax.ab2A_X,['\DeltaLFP' char(10) '(mV)'])
set(ax.ab2A_X,'XLim',[0 70],'XTick',[])
set(ax.ab2A_X,'YDir','reverse')
ylabel(ax.ab2A_R,['Firing' char(10) 'rate (Hz)'])
set(ax.ab2A_R,'XLim',[0 70])
xlabel(ax.ab2A_R,'Time (s)')


% move all the time traces to the left, and seperate the ab2A and ab3A responses 
spacing = .115;
ax.ab3A_S.Position = [.1 .82 .37 .1];
ax.ab3A_X.Position = [.1 .82 - spacing .37 .1];
ax.ab3A_R.Position = [.1 .82 - 2*spacing .37 .1];
ax.ab2A_S.Position = [.1 .32 .37 .1];
ax.ab2A_X.Position = [.1 .32 - spacing .37 .1];
ax.ab2A_R.Position = [.1 .32 - 2*spacing .37 .1];

set(ax.ab3A_SZi2,'YLim',[0 2.5])
set(ax.ab3A_SZi,'YLim',[0 2.5])
set(ax.ab3A_SZ,'YLim',[0 2.5])

figure(fig1)
prettyFig('fs',14,'plw',1.5);


% add some text saying which dataset is which
canvas = axes;
canvas.Position = [0 0 1 1];
canvas.XTick = []; 
canvas.YTick = [];
uistack(canvas,'bottom');

th = text(.1,.1,'ab3A responses to ethyl acetate');
th.Position = [.2 .95];
th.FontSize = 18;
th.Parent = canvas;


th = text(.1,.1,'ab2A responses to 2-butanone');
th.Position  = [.2 .47];
th.FontSize = 18;
th.Parent = canvas;

% add rectangles to indicate which region of the dose response we compute deviations in 
axes(ax.ab2A_drX)
h = rectangle('Position',[.3,-30,2.7,30]);
h.EdgeColor = [.95 .95 .95];
h.FaceColor = [.95 .95 .95];
uistack(h,'bottom')

pause(1)
axes(ax.ab2A_drR)
pause(1)
h = rectangle('Position',[.3,0,2.7,300]);
h.EdgeColor = [.95 .95 .95];
h.FaceColor = [.95 .95 .95];
uistack(h,'bottom')

% deintersect the dose-response axes
ax.ab2A_drX.XLim(2) = 11;
ax.ab2A_drX.YLim(2) = 0;
ax.ab2A_drR.XLim(2) = 11;
ax.ab2A_drR.YLim(1) = 0;

deintersectAxes(ax.ab2A_drX)
deintersectAxes(ax.ab2A_drR)

% label things
labelAxes(ax.ab3A_S,'a','x_offset',-.01,'font_size',24);
labelAxes(ax.ab2A_S,'b','x_offset',-.01,'font_size',24);

labelAxes(ax.X_proj,'c','x_offset',-.01,'font_size',24);
labelAxes(ax.R_proj,'d','x_offset',-.01,'font_size',24);
labelAxes(ax.ab2A_drX,'e','x_offset',-.01,'font_size',24);
labelAxes(ax.ab2A_drR,'f','x_offset',-.01,'font_size',24);

% now do fig 2
axis(ax.deviations_R_300ms,'square')
axis(ax.deviations_X_300ms,'square')
axis(ax.deviations_prev_whiff_X,'square')
axis(ax.deviations_prev_whiff_X,'square')

% whiff plots
xlabel(ax.ab2A_RZ,'Time since whiff (ms)')

% clean up the zoom plots
ax.ab3A_SZ.XTick = [];
ax.ab3A_XZ.XTick = [];
ax.ab2A_SZ.XTick = [];
ax.ab2A_XZ.XTick = [];

ax.ab3A_SZ.XLim = [-300 300];
ax.ab3A_XZ.XLim = [-300 300];
ax.ab2A_RZ.XLim = [-300 300];
ax.ab2A_SZ.XLim = [-300 300];
ax.ab2A_XZ.XLim = [-300 300];
ax.ab3A_RZ.XLim = [-300 300];

% clean up the insets a little
ax.ab3A_SZi.XTick = [];
ax.ab3A_SZi2.XTick = [];
ax.ab3A_XZi.XTick = [];
ax.ab3A_RZi.XTick = [];
ax.ab2A_SZi.XTick = [];
ax.ab2A_SZi2.XTick = [];
ax.ab2A_XZi.XTick = [];
ax.ab2A_RZi.XTick = [];
ax.ab3A_SZi.YLim = [0 2.5];
ax.ab3A_SZi.YTick = [0 1 2];
ax.ab3A_SZi2.YLim = [0 2.5];
ax.ab3A_SZi2.YTick = [0 1 2];
ax.ab3A_XZi.YLim = [-15 0];
ax.ab3A_XZi.YTick = [-15 0];
ax.ab3A_XZ.YDir  = 'reverse'; 
ax.ab3A_XZi.YDir  = 'reverse'; 
ax.ab3A_RZi.YLim = [0 160];
ax.ab3A_RZi.YTick = [0 150];
ax.ab2A_SZi.YLim = [0 1.5];
ax.ab2A_SZi2.YLim = [0 1.5];
ax.ab2A_SZi.YTick = [0 1];
ax.ab2A_SZi2.YTick = [0 1];
ax.ab2A_XZi.YLim = [-20 0];
ax.ab2A_XZi.YTick = [-20 0];
ax.ab2A_XZ.YDir  = 'reverse';
ax.ab2A_XZi.YDir  = 'reverse';  
ax.ab2A_RZi.YLim = [0 250];
ax.ab2A_RZi.YTick = [0 250];

% fix xlims for all the bar charts
fn = fieldnames(ax);
for i = 1:length(fn)
	if any(strfind(fn{i},'Zi'))
		ax.(fn{i}).XLim = [0 5];
	end
end



spacing = .115;
inset_height = .07;

ax.ab3A_SZ.Position = [.1 .79 .1 .1];
ax.ab3A_SZi.Position = [.25 .79 .05 inset_height];
ax.ab3A_SZi2.Position = [.35 .79 .05 inset_height];
ax.ab3A_XZ.Position = [.1 .79 - spacing .1 .1];
ax.ab3A_XZi.Position = [.35 .79 - spacing .05 inset_height];
ax.ab3A_RZ.Position = [.1 .79 - 2*spacing .1 .1];
ax.ab3A_RZi.Position = [.35 .79 - 2*spacing .05 inset_height];

ax.ab2A_SZ.Position = [.1 .32 .1 .1];
ax.ab2A_SZi.Position = [.25 .32 .05 inset_height];
ax.ab2A_SZi2.Position = [.35 .32 .05 inset_height];
ax.ab2A_XZ.Position = [.1 .32 - spacing .1 .1];
ax.ab2A_XZi.Position = [.35 .32 - spacing .05 inset_height];
ax.ab2A_RZ.Position = [.1 .32 - 2*spacing .1 .1];
ax.ab2A_RZi.Position = [.35 .32 - 2*spacing .05 .08];

figure(fig2)
prettyFig('fs',14,'plw',1.5);

% clean up the deviations vs. stimulus plots

ax.deviations_X_300ms.Position = [.48 .55 .2 .37];
ax.deviations_X_300ms.YLim = [-.4 .8];
lh.x_300.Location = 'southeast';
axis(ax.deviations_X_0ms,'square')
ax.deviations_X_0ms.Position = [.54 .75 .15 .15];
ax.deviations_X_0ms.Box = 'off';
ax.deviations_X_0ms.XScale = 'linear';
ax.deviations_X_0ms.XLim = [.3 3];
ax.deviations_X_0ms.YLim = [-.25 .25];
lh.x_no_corr.Position = [.55 .91 .1 .02];

ax.deviations_R_300ms.Position = [.75 .55 .2 .37];
ax.deviations_R_300ms.YLim = [-.4 .8];
lh.r_300.Location = 'southeast';
axis(ax.deviations_R_0ms,'square')

ax.deviations_R_0ms.Box = 'off';
ax.deviations_R_0ms.XScale = 'linear';
ax.deviations_R_0ms.XLim = [.3 3];
ax.deviations_R_0ms.YLim = [-.25 .25];
ax.deviations_R_0ms.Position = [.84 .75 .15 .15];
lh.r_no_corr.Position = [.865 .91 .1 .02];

% fix the T_before vs. S_before plots
ax.deviations_prev_whiff_X.Position = [.48 .075 .2 .37];
ax.deviations_prev_whiff_R.Position = [.75 .075 .2 .37];

ylabel(ax.ab3A_SZ,['Stimulus' char(10) '(V)'])
ylabel(ax.ab3A_XZ,['\DeltaLFP' char(10) '(mV)'])
ylabel(ax.ab3A_RZ,['Firing' char(10) 'rate (Hz)'])

ylabel(ax.ab2A_SZ,['Stimulus' char(10) '(V)'])
ylabel(ax.ab2A_XZ,['\DeltaLFP' char(10) '(mV)'])
ylabel(ax.ab2A_RZ,['Firing' char(10) 'rate (Hz)'])

% add some text saying which dataset is which
canvas = axes;
canvas.Position = [0 0 1 1];
canvas.XTick = []; 
canvas.YTick = [];
uistack(canvas,'bottom');

th = text(.1,.1,['ab3A responses' char(10) 'to ethyl acetate']);
th.Position = [.1 .95];
th.FontSize = 18;
th.Parent = canvas;


th = text(.1,.1,['ab2A responses' char(10) 'to 2-butanone']);
th.Position  = [.1 .47];
th.FontSize = 18;
th.Parent = canvas;

% label things
labelAxes(ax.ab3A_SZ,'a','x_offset',-.05,'font_size',24);
labelAxes(ax.ab2A_SZ,'b','x_offset',-.05,'font_size',24);

labelAxes(ax.deviations_X_300ms,'c','x_offset',-.01,'font_size',24);
labelAxes(ax.deviations_R_300ms,'d','x_offset',-.01,'font_size',24);

labelAxes(ax.deviations_prev_whiff_X,'e','x_offset',-.01,'font_size',24);
labelAxes(ax.deviations_prev_whiff_R,'f','x_offset',-.01,'font_size',24);


if being_published
	figure(fig1)
	snapnow
	delete(fig1)
	snapnow
	delete(fig2)
end

return

 ;;;;;;  ;;     ;; ;;;;;;;;  ;;;;;;;;     ;;;;;;;; ;;;;  ;;;;;;   
;;    ;; ;;     ;; ;;     ;; ;;     ;;    ;;        ;;  ;;    ;;  
;;       ;;     ;; ;;     ;; ;;     ;;    ;;        ;;  ;;        
 ;;;;;;  ;;     ;; ;;;;;;;;  ;;;;;;;;     ;;;;;;    ;;  ;;   ;;;; 
      ;; ;;     ;; ;;        ;;           ;;        ;;  ;;    ;;  
;;    ;; ;;     ;; ;;        ;;           ;;        ;;  ;;    ;;  
 ;;;;;;   ;;;;;;;  ;;        ;;           ;;       ;;;;  ;;;;;;   


 ;;;;;;  ;;;;;;;; ;;;; ;;     ;; ;;     ;; ;;       ;;     ;;  ;;;;;;  
;;    ;;    ;;     ;;  ;;;   ;;; ;;     ;; ;;       ;;     ;; ;;    ;; 
;;          ;;     ;;  ;;;; ;;;; ;;     ;; ;;       ;;     ;; ;;       
 ;;;;;;     ;;     ;;  ;; ;;; ;; ;;     ;; ;;       ;;     ;;  ;;;;;;  
      ;;    ;;     ;;  ;;     ;; ;;     ;; ;;       ;;     ;;       ;; 
;;    ;;    ;;     ;;  ;;     ;; ;;     ;; ;;       ;;     ;; ;;    ;; 
 ;;;;;;     ;;    ;;;; ;;     ;;  ;;;;;;;  ;;;;;;;;  ;;;;;;;   ;;;;;;  

 ;;;;;;  ;;;;;;;;    ;;;    ;;;;;;;; ;;;;  ;;;;;;  ;;;;;;;; ;;;;  ;;;;;;   ;;;;;;  
;;    ;;    ;;      ;; ;;      ;;     ;;  ;;    ;;    ;;     ;;  ;;    ;; ;;    ;; 
;;          ;;     ;;   ;;     ;;     ;;  ;;          ;;     ;;  ;;       ;;       
 ;;;;;;     ;;    ;;     ;;    ;;     ;;   ;;;;;;     ;;     ;;  ;;        ;;;;;;  
      ;;    ;;    ;;;;;;;;;    ;;     ;;        ;;    ;;     ;;  ;;             ;; 
;;    ;;    ;;    ;;     ;;    ;;     ;;  ;;    ;;    ;;     ;;  ;;    ;; ;;    ;; 
 ;;;;;;     ;;    ;;     ;;    ;;    ;;;;  ;;;;;;     ;;    ;;;;  ;;;;;;   ;;;;;;  

cdata = consolidateData2(getPath(dataManager,'4608c42b12191b383c84fea52392ea97'));


figure('outerposition',[0 0 1200 801],'PaperUnits','points','PaperSize',[1200 801]); hold on

subplot(2,3,1); hold on

% plot pdf of whiff intensities 

% show whiff statistics 
S = cdata.PID(:,cdata.orn==4); 
S = bsxfun(@minus,S,min(S));
S = S(:);
ws = whiffStatistics(S,0*S,0*S,300,'MinPeakProminence',max(S/1e2),'debug',false);


[y,x] = histcounts(vertcat(ws.stim_peaks),30);
y = y/sum(y); x = x(2:end)+mean(diff(x));
plot(x,y,'k+')

m = fittype('log((a./x).*exp(-x./b))');
ff = fit(x(y>0)',log(y(y>0))',m,'Upper',[10 100],'Lower',[0 1],'StartPoint',[.14 10.37]);
x = logspace(-1,2,100);
plot(sort(x),exp(ff(sort(x))),'r')

set(gca,'XScale','log','YScale','log','XTick',[.1 1 10 100],'YTick',[1e-3 1e-2 1e-1 1],'XLim',[.1 100],'YLim',[1e-3 1])
xlabel('Whiff intensity (V)')
ylabel('Probability')

th(2) = text(5, .5,'$\sim\frac{1}{c}\exp\left(-\frac{c}{C}\right)$','interpreter','latex','Color','r','FontSize',20);

% plot whiff durations 
[ons,offs] = computeOnsOffs(S>.01);
wd = offs - ons; wd(wd<10) = []; % whiffs this brief are artifacts 

subplot(2,3,2); hold on
[y,x] = histcounts(wd,50);
y = y/sum(y); x = x(2:end)+mean(diff(x));
plot(x,y,'k+')
set(gca,'XScale','log','YScale','log')

a = 1; m = fittype('a + n*x');
xx = vectorise(log(x)); yy = vectorise(log(y));
ff = fit(xx(yy>-Inf),yy(yy>-Inf),m,'Upper',[Inf -1.5],'Lower',[-Inf -1.5],'StartPoint',[6 -1.5]);
plot(x,exp(ff(log(x))),'r')
xlabel('Whiff duration (ms)')
ylabel('Probability')
th(3) = text(1e3, .1,'$\sim t_{w}^{-\frac{3}{2}}$','interpreter','latex','Color','r','FontSize',20);

% plot blank durations 
[ons,offs] = computeOnsOffs(S<.01);
wd = offs - ons; wd(wd<10) = []; % whiffs this brief are artifacts 

subplot(2,3,3); hold on
[y,x] = histcounts(wd,50);
y = y/sum(y); x = x(2:end)+mean(diff(x));
plot(x,y,'k+')
set(gca,'XScale','log','YScale','log','YLim',[1e-3 1])

a = 1; m = fittype('a + n*x');
xx = vectorise(log(x)); yy = vectorise(log(y));
ff = fit(xx(yy>-Inf),yy(yy>-Inf),m,'Upper',[Inf -1.5],'Lower',[-Inf -1.5],'StartPoint',[6 -1.5]);
plot(x,exp(ff(log(x))),'r')
xlabel('Blank duration (ms)')
ylabel('Probability')
th(3) = text(1e4, .1,'$\sim t_{b}^{-\frac{3}{2}}$','interpreter','latex','Color','r','FontSize',20);

% now compute the correlation between the mean and variance across many windows 
window_sizes = factor2(length(S));
window_sizes = window_sizes(window_sizes < 1e4 & window_sizes > 10);
r2 = NaN*window_sizes;
for i = 1:length(window_sizes)
	tS = reshape(S,window_sizes(i),length(S)/window_sizes(i));
	r2(i) = rsquare(mean(tS),std(tS));
end

subplot(2,3,5); hold on
plot(window_sizes,r2,'k+')
xlabel('Window size (ms)')
ylabel('r^2_{\mu,\sigma}')
set(gca,'XScale','log','XTick',[10 100 1e3 1e4])

subplot(2,3,4); hold on
ws = window_sizes(27);
tS = reshape(S,ws,length(S)/ws);
plot(mean(tS),std(tS),'.','Color',[.3 .3 .3])
plot([1e-4 100],[1e-4 100],'k--')
set(gca,'XScale','log','YScale','log','XLim',[1e-2 10],'YLim',[1e-2 10],'XTick',[1e-2 1e-1 1 10])
xlabel('\mu_{S} in preceding 400ms(V)')
ylabel('\sigma_{S} in preceding 400ms (V)')

% now show the autocorrelation function
subplot(2,3,6); hold on
S = cdata.PID(:,cdata.orn==4); 
S = bsxfun(@minus,S,min(S));
for i = 1:size(S,2)
	[ac(:,i),lags] = autocorr(S(:,i),length(S)-1);
end
errorShade(lags,mean(ac,2),std(ac,[],2),'Color',[0 0 0]);
set(gca,'XScale','log','XTick',[1 10 100 1e3 1e4 1e5])
xlabel('Lag (ms)')
ylabel('Autocorrelation')

prettyFig();
labelFigure('x_offset',0)

if being_published
	snapnow
	delete(gcf)
end


 ;;;;;;  ;;     ;; ;;;;;;;;  ;;;;;;;;     ;;;;;;;; ;;;;  ;;;;;;   
;;    ;; ;;     ;; ;;     ;; ;;     ;;    ;;        ;;  ;;    ;;  
;;       ;;     ;; ;;     ;; ;;     ;;    ;;        ;;  ;;        
 ;;;;;;  ;;     ;; ;;;;;;;;  ;;;;;;;;     ;;;;;;    ;;  ;;   ;;;; 
      ;; ;;     ;; ;;        ;;           ;;        ;;  ;;    ;;  
;;    ;; ;;     ;; ;;        ;;           ;;        ;;  ;;    ;;  
 ;;;;;;   ;;;;;;;  ;;        ;;           ;;       ;;;;  ;;;;;;   

  ;;;;;;      ;;;    ;;     ;;  ;;;;;;   ;;;;;;  ;;;;    ;;;    ;;    ;; 
;;    ;;    ;; ;;   ;;     ;; ;;    ;; ;;    ;;  ;;    ;; ;;   ;;;   ;; 
;;         ;;   ;;  ;;     ;; ;;       ;;        ;;   ;;   ;;  ;;;;  ;; 
;;   ;;;; ;;     ;; ;;     ;;  ;;;;;;   ;;;;;;   ;;  ;;     ;; ;; ;; ;; 
;;    ;;  ;;;;;;;;; ;;     ;;       ;;       ;;  ;;  ;;;;;;;;; ;;  ;;;; 
;;    ;;  ;;     ;; ;;     ;; ;;    ;; ;;    ;;  ;;  ;;     ;; ;;   ;;; 
 ;;;;;;   ;;     ;;  ;;;;;;;   ;;;;;;   ;;;;;;  ;;;; ;;     ;; ;;    ;; 

 ;;;;;;;; ;;;; ;;       ;;;;;;;; ;;;;;;;; ;;;;;;;;   ;;;;;;  
;;        ;;  ;;          ;;    ;;       ;;     ;; ;;    ;; 
;;        ;;  ;;          ;;    ;;       ;;     ;; ;;       
;;;;;;    ;;  ;;          ;;    ;;;;;;   ;;;;;;;;   ;;;;;;  
;;        ;;  ;;          ;;    ;;       ;;   ;;         ;; 
;;        ;;  ;;          ;;    ;;       ;;    ;;  ;;    ;; 
;;       ;;;; ;;;;;;;;    ;;    ;;;;;;;; ;;     ;;  ;;;;;;  



cdata = consolidateData2(getPath(dataManager,'4608c42b12191b383c84fea52392ea97'));
[cdata, data] =  assembleScaledNatStim(cdata);
time = 1e-3*(1:length(data(1).S));

% compute filters from this neuron
clear K1 K2
for j = 1:size(data(2).S,2)
	S = data(2).S(:,j);
	X = data(2).X(:,j);
	R = data(2).R(:,j);

	temp1 = fitFilter2Data(S,X,'offset',200);
	temp2 = fitFilter2Data(S,R,'offset',200);
	K1(:,j) = temp1(100:end-100);
	K2(:,j) = temp2(100:end-100);
end
filtertime = 1e-3*(1:length(K1)) - .1;

i = 2;
clear Xp Rp
for j = 1:size(data(i).S,2)
	S = data(i).S(:,j);
	X = data(i).X(:,j);
	R = data(i).R(:,j);

	Xp(:,j) = convolve(time,S,K1(:,j),filtertime);
	Rp(:,j) = convolve(time,S,K2(:,j),filtertime);
	K1(:,j) = K1(:,j)/(nanstd(Xp(:,j))/nanstd(S)); % normalise correctly 
	Xp(:,j) = convolve(time,S,K1(:,j),filtertime);
	K2(:,j) = K2(:,j)/(nanstd(Rp(:,j))/nanstd(S)); % normalise correctly 
	Rp(:,j) = convolve(time,S,K2(:,j),filtertime);
end
nat.Rp = Rp;
nat.Xp = Xp;

% now get filters from the gaussian data
MSGdata = consolidateData2(getPath(dataManager,'3ea08ccfa892c6545d74bbdaaa6cbee1'));
MSGdata = cleanMSGdata(MSGdata,'extract_filter',true);

gK1 = MSGdata.K1(:,MSGdata.paradigm<3);
gK2 = MSGdata.K2(:,MSGdata.paradigm<3);

nK1 = nanmean(gK1,2);
nK2 = nanmean(gK2,2);

i = 2;
clear Xp Rp
filtertime = (1:length(nK1))*1e-3 - .1;
for j = 1:size(data(i).S,2)
	S = data(i).S(:,j);
	X = data(i).X(:,j);
	R = data(i).R(:,j);

	Xp(:,j) = convolve(time,S,nK1,filtertime);
	Rp(:,j) = convolve(time,S,nK2,filtertime);
	nK1 = nK1/(nanstd(Xp(:,j))/nanstd(S)); % normalise correctly 
	Xp(:,j) = convolve(time,S,nK1,filtertime);
	nK2 = nK2/(nanstd(Rp(:,j))/nanstd(S)); % normalise correctly 
	Rp(:,j) = convolve(time,S,nK2,filtertime);

end



figure('outerposition',[0 0 1300 901],'PaperUnits','points','PaperSize',[1300 901]); hold on
clear ax
for i = 1:6
	ax(i) = subplot(2,3,i); hold on
end
% plot the filers and the projections 
K1 = K1/norm(K1);
filtertime = (1:length(K1))*1e-3 - .1;
errorShade(ax(1),filtertime,nanmean(K1,2),sem(K1'),'Color','k');

K2 = K2/norm(K2);
errorShade(ax(4),filtertime,nanmean(K2,2),sem(K2'),'Color','k');

gK1 = gK1/norm(gK1);
filtertime = (1:length(gK1))*1e-3 - .1;
errorShade(ax(1),filtertime,nanmean(gK1,2),sem(gK1'),'Color','r');

gK2 = gK2/norm(gK2);
errorShade(ax(4),filtertime,nanmean(gK2,2),sem(gK2'),'Color','r');


% now plot the projections
l = plot(ax(2),nat.Xp(:),data(2).X(:),'Color','k');
r2 = rsquare(nat.Xp(:),data(2).X(:));
legend(l,['r^2 = ' oval(r2)],'Location','southeast');
l = plot(ax(5),nat.Rp(:),data(2).R(:),'Color','k');
r2 = rsquare(nat.Rp(:),data(2).R(:));
legend(l,['r^2 = ' oval(r2)],'Location','southeast');

l = plot(ax(3),Xp(:),data(2).X(:),'Color','r');
r2 = rsquare(Xp(:),data(2).X(:));
legend(l,['r^2 = ' oval(r2)],'Location','southeast');
l = plot(ax(6),Rp(:),data(2).R(:),'Color','r');
r2 = rsquare(Rp(:),data(2).X(:));
legend(l,['r^2 = ' oval(r2)],'Location','southeast');

% label axes, etc
xlabel(ax(1),'Filter lag (s)')
ylabel(ax(1),'LFP filter (norm)')
set(ax(1:3),'YDir','reverse')
set(ax(2:3),'XDir','reverse')

xlabel(ax(4),'Filter lag (s)')
ylabel(ax(4),'Firing filter (norm)')

xlabel(ax(2),'Projected stimulus (V)')
xlabel(ax(3),'Projected stimulus (V)')
xlabel(ax(5),'Projected stimulus (V)')
xlabel(ax(6),'Projected stimulus (V)')

ylabel(ax(2),'\DeltaLFP (mV)')
ylabel(ax(3),'\DeltaLFP (mV)')
ylabel(ax(5),'Firing rate (Hz)')
ylabel(ax(6),'Firing rate (Hz)')

prettyFig();

labelFigure('x_offset',0)

shrinkDataInPlot(ax([2 3 5 6]),5)

if being_published
	snapnow
	delete(gcf)
end


clearvars -except being_published

;;    ;; ;;          ;;     ;;  ;;;;;;;  ;;;;;;;;  ;;;;;;;; ;;       
;;;   ;; ;;          ;;;   ;;; ;;     ;; ;;     ;; ;;       ;;       
;;;;  ;; ;;          ;;;; ;;;; ;;     ;; ;;     ;; ;;       ;;       
;; ;; ;; ;;          ;; ;;; ;; ;;     ;; ;;     ;; ;;;;;;   ;;       
;;  ;;;; ;;          ;;     ;; ;;     ;; ;;     ;; ;;       ;;       
;;   ;;; ;;          ;;     ;; ;;     ;; ;;     ;; ;;       ;;       
;;    ;; ;;;;;;;;    ;;     ;;  ;;;;;;;  ;;;;;;;;  ;;;;;;;; ;;;;;;;; 


% supp figure
% Can a NL model fit to the LFP show context-dependent response modulation? 

cdata = consolidateData2(getPath(dataManager,'4608c42b12191b383c84fea52392ea97'));
[cdata, data] =  assembleScaledNatStim(cdata);

clear fd
for i = 1:3
	fd(i).response = data(2).X(:,i);
	S = data(2).S(:,i); S = S - min(S);
	fd(i).stimulus = [S fd(i).response];
end



clear p
p.k_D = 0.20742;
p.n = 1.475;

% generate predictions
clear K
for i = 1:3
	S = data(2).S(:,i); S = S - min(S);
	R = data(2).X(:,i); R = R - min(R);
	[data(2).XP(:,i) ,K(:,i)] = NLModel([S R],p);
	K(:,i) = K(:,i)/norm(K(:,i));
end


figure('outerposition',[0 0 1403 801],'PaperUnits','points','PaperSize',[1403 801]); hold on
% first show the best fit model
x = logspace(-2,2,100);
y = 1./(1+(p.k_D./x).^p.n);
ax(1) = subplot(2,4,1); hold on
plot(x,y,'r')
set(gca,'XScale','log')
xlabel('Stimulus')
ylabel('a')


ax(2) = subplot(2,4,2); hold on
filtertime = 1:length(K); filtertime = filtertime - 50;
errorShade(filtertime,mean(K,2),std(K,[],2),'Color','r');
xlabel('Filter lag (ms)')
ylabel('Filter (norm)')
set(gca,'YDir','reverse')

% show one trace
example_trace = 3;
ax(3) = subplot(2,4,3:4); hold on
time = 1:length(data(2).S); time = time*1e-3;
clear l
l(1) = plot(time,data(2).X(:,example_trace),'k');
l(2) = plot(time,data(2).XP(:,example_trace),'r');
set(gca,'XLim',[0 70],'YDir','reverse')
xlabel('Time (s)')
ylabel('\Delta LFP (mV)')
r2 = rsquare(data(2).X(:,example_trace),data(2).XP(:,example_trace));
legend(l,{'ab2 LFP',['NL model, r^2 = ' oval(r2)]},'Location','northwest')


% now pull out a filter from this
Khat = NaN*K;
for i = 1:3
	S = data(2).S(:,i); S = S - min(S);
	R = data(2).X(:,i); R = R - min(R);
	temp = fitFilter2Data(S,R,'filter_length',700,'offset',100);
	Khat(:,i) = temp(50:end-50);
	Khat(:,i) = Khat(:,i)/norm(Khat(:,i));
end

ax(5) = subplot(2,4,5); hold on
filtertime = 1:length(K); filtertime = filtertime - 50;
errorShade(filtertime,mean(Khat,2),std(Khat,[],2),'Color','k');
xlabel('Filter lag (ms)')
ylabel('Filter K (norm)')
set(gca,'YDir','reverse')

% make projections using this 
for i = 1:3
	data(2).fp(:,i) = (1/3)*convolve(1:length(data(2).S(:,i)),data(2).S(:,i),Khat(:,i),filtertime);
end

ax(6) = subplot(2,4,6); hold on
l = plot(data(2).fp(:),data(2).XP(:),'k');
r2 = ['r^2 = ' oval(rsquare(data(2).fp(:),data(2).XP(:)))];
legend(l,r2,'Location','southeast')
xlabel('K \otimes s(t)')
ylabel('NL model prediction (\Delta mV)')
set(gca,'XDir','reverse','YDir','reverse','XLim',[-20 0 ])


ax(7) = subplot(2,4,7); hold on
ax(8) = subplot(2,4,8); hold on
show_these = [2       27764
           2       28790
           2       59776
           3       58049];

for i = 2
	for j = 1:length(show_these)
		this_stim = show_these(j,1);
		this_loc = show_these(j,2);

		S = data(i).S(:,this_stim);
		X = data(i).X(:,this_stim);
		XP = data(i).XP(:,this_stim);

		a = this_loc - 300;
		z = this_loc+300;

		t = (1:length(S(a:z))) - 300;
		plot(ax(7),t,X(a:z))
		plot(ax(8),t,XP(a:z))

	end
end

xlabel(ax(7),'Time since whiff (ms)')
ylabel(ax(7),'ab2 LFP (mV)')
set(ax(7),'YDir','reverse')

xlabel(ax(8),'Time since whiff (ms)')
ylabel(ax(8),'NL model (mV)')
set(ax(8),'YDir','reverse')

prettyFig();
labelFigure('x_offset',0)


if being_published
	snapnow
	delete(gcf)
end



%% Version Info
%
pFooter;


