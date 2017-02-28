% fig1_nat_stim_example
% now uses ab2A and ab3A data 

%%
% This script makes fig 1 for the eLife resubmission. The point of this figure is to show our naturalistic stimulus data, and the features of the ORN response we see in it. 

pHeader;


% make the figure and all subplots 
figure('outerposition',[0 0 900 1001],'PaperUnits','points','PaperSize',[900 1001]); hold on
clear ax

% time series for ab3A
ax.ab3A_S = subplot(7,4,1:3); hold on
ax.ab3A_X = subplot(7,4,5:7); hold on
ax.ab3A_R = subplot(7,4,9:11); hold on

% time series for ab2A
ax.ab2A_S = subplot(7,4,13:15); hold on
ax.ab2A_X = subplot(7,4,17:19); hold on
ax.ab2A_R = subplot(7,4,21:23); hold on

% zooms into time series, one for ab3A, one for ab2A 
ax.ab3A_SZ = subplot(7,4,4); hold on
ax.ab3A_XZ = subplot(7,4,8); hold on
ax.ab3A_RZ = subplot(7,4,12); hold on

ax.ab2A_SZ = subplot(7,4,16); hold on
ax.ab2A_XZ = subplot(7,4,20); hold on
ax.ab2A_RZ = subplot(7,4,24); hold on

% filters + projections for ab2A 
ax.X_proj = subplot(7,4,25); hold on
ax.R_proj = subplot(7,4,26); hold on

% dose response for ab2A
ax.ab2A_drX = subplot(7,4,27); hold on
ax.ab2A_drR = subplot(7,4,28); hold on

% add insets showing that whiffs don't change, but responses do
ax.ab3A_SZi = axes; hold on
ax.ab3A_SZi.Position = [.88 .89 .05 .05];

ax.ab3A_XZi = axes; hold on
ax.ab3A_XZi.Position = [.88 .75 .05 .05];

ax.ab3A_RZi = axes; hold on
ax.ab3A_RZi.Position = [.88 .63 .05 .05];

ax.ab2A_SZi = axes; hold on
ax.ab2A_SZi.Position = [.88 .51 .05 .05];

ax.ab2A_XZi = axes; hold on
ax.ab2A_XZi.Position = [.88 .39 .05 .05];

ax.ab2A_RZi = axes; hold on
ax.ab2A_RZi.Position = [.88 .28 .05 .05];


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

for i = 1:length(show_these)
	this_paradigm = show_these(i,1); 

	PID = data(2).S(:,this_paradigm);
	LFP = data(2).X(:,this_paradigm);
	fA =  data(2).R(:,this_paradigm);

	temp = nanmean(PID(show_these(i,2):show_these(i,2)+10,:));
	S_int(i) = nanmean(temp);
	S_int_err(i) = nanstd(temp)/sqrt(length(temp));

	temp = nanmean(LFP(show_these(i,2)+90:show_these(i,2)+100,:));
	X_int(i) = nanmean(temp);
	X_int_err(i) = nanstd(temp)/sqrt(length(temp));

	temp = nanmean(fA(show_these(i,2)+60:show_these(i,2)+80,:));
	R_int(i) = nanmean(temp);
	R_int_err(i) = nanstd(temp)/sqrt(length(temp));

end

cla(ax.ab2A_SZi)
cla(ax.ab2A_XZi)
cla(ax.ab2A_RZi)
c = lines(10);

[~,~,S_order] = unique(S_int);
[~,~,X_order] = unique(X_int); X_order = 5 - X_order;
[~,~,R_order] = unique(R_int);

for i = 1:length(show_these)
	superbar(ax.ab2A_SZi,S_order(i),S_int(i),'E',S_int_err(i),'BarFaceColor',c(i,:),'ErrorbarColor',c(i,:))
	superbar(ax.ab2A_XZi,X_order(i),X_int(i),'E',X_int_err(i),'BarFaceColor',c(i,:),'ErrorbarColor',c(i,:))
	superbar(ax.ab2A_RZi,R_order(i),R_int(i),'E',R_int_err(i),'BarFaceColor',c(i,:),'ErrorbarColor',c(i,:))

end

% show response vs. projected stimulus
K = fitFilter2Data(data(2).S(:,1),data(2).X(:,1),'offset',200); 
K = K(100:end-100);
filtertime = (1:length(K)) - 100;
fp = convolve(1:length(data(2).S(:)),data(2).S(:),K,filtertime);
plot(ax.X_proj,fp,data(2).X(:),'k')
set(ax.X_proj,'YDir','reverse','XDir','reverse','XLim',[-4 0],'YLim',[-30 0])
xlabel(ax.X_proj,'Projected stimulus (V)')
ylabel(ax.X_proj,'ab2 LFP (mV)')

K = fitFilter2Data(data(2).S(:,1),data(2).R(:,1),'offset',200); 
K = K(100:end-100);
filtertime = (1:length(K)) - 100;
fp = convolve(1:length(data(2).S(:)),data(2).S(:),K,filtertime);
plot(ax.R_proj,fp,data(2).R(:),'k')
set(ax.R_proj,'XLim',[0 5],'YLim',[0 300])
xlabel(ax.R_proj,'Projected stimulus (V)')
ylabel(ax.R_proj,'ab2A firing rate (Hz)')

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
show_these = [6112 26669 28413 64360];

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

for i = 1:length(show_these)
	temp = nanmean(PID(show_these(i):show_these(i)+10,:));
	S_int(i) = nanmean(temp);
	S_int_err(i) = nanstd(temp)/sqrt(length(temp));

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

[~,~,S_order] = unique(S_int);
[~,~,X_order] = unique(X_int); X_order = max(X_order)+1 - X_order;
[~,~,R_order] = unique(R_int);

for i = 1:length(show_these)
	superbar(ax.ab3A_SZi,S_order(i),S_int(i),'E',S_int_err(i),'BarFaceColor',c(i,:),'ErrorbarColor',c(i,:))
	superbar(ax.ab3A_XZi,X_order(i),X_int(i),'E',X_int_err(i),'BarFaceColor',c(i,:),'ErrorbarColor',c(i,:))
	superbar(ax.ab3A_RZi,R_order(i),R_int(i),'E',R_int_err(i),'BarFaceColor',c(i,:),'ErrorbarColor',c(i,:))

end


prettyFig('fs',12,'plw',1.5);

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

set(ax.ab3A_XZ,'YDir','reverse')
set(ax.ab2A_XZ,'YDir','reverse')
xlabel(ax.ab2A_RZ,'Time since whiff (ms)')


xlabel(ax.ab2A_drX,'Whiff amplitude (V)')
ylabel(ax.ab2A_drX,'ab2 LFP (mV)')
set(ax.ab2A_drX,'XScale','log','XLim',[1e-2 1e1],'YDir','reverse','XTick',[1e-2 1e-1 1e0 1e1 1e2])

xlabel(ax.ab2A_drR,'Whiff amplitude (V)')
ylabel(ax.ab2A_drR,'ab2A firing rate (Hz)')
set(ax.ab2A_drR,'XScale','log','XLim',[1e-2 1e1],'XTick',[1e-2 1e-1 1e0 1e1 1e2])

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
ax.ab3A_XZi.XTick = [];
ax.ab3A_RZi.XTick = [];
ax.ab2A_SZi.XTick = [];
ax.ab2A_XZi.XTick = [];
ax.ab2A_RZi.XTick = [];

ax.ab3A_SZi.YLim = [0 0.4];
ax.ab3A_SZi.YTick = [0 0.4];

ax.ab3A_XZi.YLim = [-15 0];
ax.ab3A_XZi.YTick = [-15 0];
ax.ab3A_XZi.YDir  = 'reverse'; 

ax.ab3A_RZi.YLim = [0 160];
ax.ab3A_RZi.YTick = [0 150];

ax.ab2A_SZi.YLim = [0 1];
ax.ab2A_SZi.YTick = [0 1];
ax.ab2A_XZi.YLim = [-20 0];
ax.ab2A_XZi.YTick = [-20 0];
ax.ab2A_XZi.YDir  = 'reverse'; 
ax.ab2A_RZi.YLim = [0 250];
ax.ab2A_RZi.YTick = [0 250];

% move things around
ax.ab3A_S.Position = [0.1300 0.88 0.57 0.07];
ax.ab3A_X.Position = [0.1300 0.80 0.57 0.07];
ax.ab3A_R.Position = [0.1300 0.72 0.57 0.07];
ax.ab2A_S.Position = [0.1300 0.56 0.57 0.07];
ax.ab2A_X.Position = [0.1300 0.48 0.57 0.07];
ax.ab2A_R.Position = [0.1300 0.4 0.57 0.07];

ax.ab3A_SZ.Position = [0.76 0.88 0.13 0.07];
ax.ab3A_XZ.Position = [0.76 0.80 0.13 0.07];
ax.ab3A_RZ.Position = [0.76 0.72 0.13 0.07];
ax.ab2A_SZ.Position = [0.76 0.56 0.13 0.07];
ax.ab2A_XZ.Position = [0.76 0.48 0.13 0.07];
ax.ab2A_RZ.Position = [0.76 0.40 0.13 0.07];

ax.ab3A_SZi.Position = [0.93 0.88 0.05 0.05];
ax.ab3A_XZi.Position = [0.93 0.80 0.05 0.05];
ax.ab3A_RZi.Position = [0.93 0.72 0.05 0.05];
ax.ab2A_SZi.Position = [0.93 0.56 0.05 0.05];
ax.ab2A_XZi.Position = [0.93 0.48 0.05 0.05];
ax.ab2A_RZi.Position = [0.93 0.40 0.05 0.05];


% make the bottom four plots bigger
ax.X_proj.Position = [0.1 0.15 .16 .16];
ax.R_proj.Position = [0.33 0.15 .16 .16];

ax.ab2A_drX.Position = [0.57 0.15 .16 .16];
ax.ab2A_drR.Position = [0.8 0.15 .16 .16];


% label things
labelAxes(ax.ab3A_S,'a','x_offset',-.05,'font_size',24);
labelAxes(ax.ab2A_S,'b','x_offset',-.05,'font_size',24);

labelAxes(ax.X_proj,'c','x_offset',-.025,'font_size',24);
labelAxes(ax.R_proj,'d','x_offset',-.025,'font_size',24);

labelAxes(ax.ab3A_SZ,'e','x_offset',-.025,'font_size',24);
labelAxes(ax.ab2A_SZ,'f','x_offset',-.025,'font_size',24);

labelAxes(ax.ab2A_drX,'g','x_offset',-.025,'font_size',24);
labelAxes(ax.ab2A_drR,'h','x_offset',-.025,'font_size',24);


% deintersect the dose-response axes

ax.ab2A_drX.XLim(2) = 11;
ax.ab2A_drX.YLim(2) = 0;
ax.ab2A_drR.XLim(2) = 11;
ax.ab2A_drR.YLim(1) = 0;

deintersectAxes(ax.ab2A_drX)
deintersectAxes(ax.ab2A_drR)

% add some text saying which dataset is which
canvas = axes;
canvas.Position = [0 0 1 1];
canvas.XTick = []; 
canvas.YTick = [];
uistack(canvas,'bottom');

th = text(.1,.1,'ab3A responses to ethyl acetate');
th.Position  = [.3 .97];
th.FontSize = 18;
th.Parent = canvas;


th = text(.1,.1,'ab2A responses to 2-butanone');
th.Position  = [.3 .65];
th.FontSize = 18;
th.Parent = canvas;

return

% shrink data in plot
shrinkDataInPlot(ax.ab2A_S,2)
shrinkDataInPlot(ax.ab2A_X,2)
shrinkDataInPlot(ax.ab2A_R,2)

shrinkDataInPlot(ax.ab3A_S,2)
shrinkDataInPlot(ax.ab3A_X,2)
shrinkDataInPlot(ax.ab3A_R,2)


shrinkDataInPlot(ax.X_proj,2)
shrinkDataInPlot(ax.R_proj,2)

if being_published
	snapnow
	delete(gcf)
end




return



















% show linear analysis output plot


ss = 1;

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
clear Xl Rl XL RL
for j = 1:size(data(i).S,2)
	S = data(i).S(:,j);
	X = data(i).X(:,j);
	R = data(i).R(:,j);

	Xp = convolve(time,S,K1(:,j),filtertime);
	Rp = convolve(time,S,K2(:,j),filtertime);
	K1(:,j) = K1(:,j)/(nanstd(Xp)/nanstd(S)); % normalise correctly 
	Xp = convolve(time,S,K1(:,j),filtertime);
	K2(:,j) = K2(:,j)/(nanstd(Rp)/nanstd(S)); % normalise correctly 
	Rp = convolve(time,S,K2(:,j),filtertime);

 
	Xl(j) = plot(ax(9),Xp(1:ss:end),X(1:ss:end),'Color',c(j,:));
	XL{j} = ['r^2 = ' oval(rsquare(Xp,X))];
	Rl(j) = plot(ax(10),Rp(1:ss:end),R(1:ss:end),'Color',c(j,:));
	RL{j} = ['r^2 = ' oval(rsquare(Rp,R))];

end
legend(Xl,XL,'Location','southeast')
legend(Rl,RL,'Location','southeast')


% axes modifications, labels, etc.
ylabel(ax(1),'Stimulus (V)')
set(ax(1),'XLim',[0 70])
ylabel(ax(2),'\Delta LFP (mV)')
set(ax(2),'XLim',[0 70])
ylabel(ax(3),['ab2A firing' char(10) 'rate (Hz)'])
xlabel(ax(3),'Time (s)')
set(ax(3),'XLim',[0 70])

xlabel(ax(4),'Whiff intensity (V)')
ylabel(ax(4),'Peak LFP (mV)')
set(ax(4),'XScale','log','XLim',[1e-2 1e1],'XTick',[1e-2 1e-1 1e0 1e1])
axes(ax(4))
axis square

xlabel(ax(5),'Whiff intensity (V)')
ylabel(ax(5),'Peak firing rate (Hz)')
set(ax(5),'XScale','log','YLim',[0 300],'XLim',[1e-2 1e1],'XTick',[1e-2 1e-1 1e0 1e1])
axes(ax(5))
axis square

xlabel(ax(8),'Time since whiff (ms)')

xlabel(ax(9),'Projected stimulus (V)')
ylabel(ax(9),'LFP (mV)')
set(ax(9),'YDir','reverse','XDir','reverse')
axes(ax(9))
axis square

xlabel(ax(10),'Projected stimulus (V)')
ylabel(ax(10),'Firing rate (Hz)')
axes(ax(10))
axis square


% move some things around
ax(1).Position = [0.1300    0.8    0.5689    0.11];
ax(2).Position = [ 0.1300    0.625    0.5689    0.11];
ax(3).Position = [ 0.1300    0.475    0.5689    0.11];
for i = 6:8
	ax(i).Position(2) = ax(i-5).Position(2);
	ax(i).Position(4) = ax(i-5).Position(4);
end

for i = [4 5 9 10]
	movePlot(ax(i),'left',.04)
	ax(i).Position(2) = .15;
	ax(i).Position(3:4) = .22;
end

% flip all the LFP axes
set(ax(4),'YDir','reverse')
set(ax(2),'YDir','reverse')
set(ax(7),'YDir','reverse')

prettyFig()

shrinkDataInPlot(ax(1),2)
shrinkDataInPlot(ax(2),2)
shrinkDataInPlot(ax(3),2)

shrinkDataInPlot(ax(9),10)
shrinkDataInPlot(ax(10),10)

% label things
labels = char(97:110);
for i = 1:length(ax)
	labelAxes(ax(i),labels(i),'x_offset',0,'y_offset',0.01,'font_size',24);
end

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


%% Supplement 1 
% In this supplement, we show the filters, and estabilish the linear analysis that we do is not complicated by filter extraction in these conditions, and that we see the same behaviour more gnerally, on other neurons. 


% now get filters from the gaussian data
MSGdata = consolidateData2(getPath(dataManager,'3ea08ccfa892c6545d74bbdaaa6cbee1'));
MSGdata = cleanMSGdata(MSGdata,'extract_filter',true);

gK1 = MSGdata.K1(:,MSGdata.paradigm<3);
gK2 = MSGdata.K2(:,MSGdata.paradigm<3);

% make the figure
figure('outerposition',[0 0 1401 902],'PaperUnits','points','PaperSize',[1401 902]); hold on
clear ax
for i = 1:12
	ax(i) = subplot(3,4,i); hold on
end

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


% normalise the filters for display
nK1 = K1; nK2 = K2;
for i = 1:3
	nK1(:,i) = K1(:,i)/norm(K1(:,i));
	nK2(:,i) = K2(:,i)/norm(K2(:,i));
end

filtertime = 1e-3*(1:length(nK1)) - .1;
errorShade(ax(1),filtertime,mean(nK1,2),sem(nK1'),'Color','k');
errorShade(ax(2),filtertime,mean(nK2,2),sem(nK2'),'Color','k');


% also show the filters, after normalising 
for i = 1:size(gK1,2)
	gK1(:,i) = gK1(:,i)/norm(gK1(:,i));
	gK2(:,i) = gK2(:,i)/norm(gK2(:,i));
end

filtertime = 1e-3*(1:length(gK1)) - .1;
errorShade(ax(1),filtertime,mean(gK1,2),sem(gK1'),'Color','r');
errorShade(ax(2),filtertime,mean(gK2,2),sem(gK2'),'Color','r');


% now project stimulus using the gaussian filters, and see how well it does
gK1m = mean(gK1,2);
gK2m = mean(gK2,2);
i = 2;
ss = 10;
clear Xl Rl XL RL
for j = 1:size(data(i).S,2)
	S = data(i).S(:,j);
	X = data(i).X(:,j);
	R = data(i).R(:,j);


	Xp = convolve(time,S,gK1m,filtertime);
	Rp = convolve(time,S,gK2m,filtertime);
	gK1m = gK1m/(nanstd(Xp)/nanstd(S)); % normalise correctly 
	Xp = convolve(time,S,gK1m,filtertime);
	gK2m = gK2m/(nanstd(Rp)/nanstd(S)); % normalise correctly 
	Rp = convolve(time,S,gK2m,filtertime);

 
	Xl(j) = plot(ax(3),Xp(1:ss:end),X(1:ss:end),'Color',c(j,:));
	XL{j} = ['r^2 = ' oval(rsquare(Xp,X))];
	Rl(j) = plot(ax(4),Rp(1:ss:end),R(1:ss:end),'Color',c(j,:));
	RL{j} = ['r^2 = ' oval(rsquare(Rp,R))];

end
legend(Xl,XL,'Location','southeast')
legend(Rl,RL,'Location','southeast')



   ;;;    ;;    ;;  ;;;;;;;  ;;;;;;;; ;;     ;; ;;;;;;;; ;;;;;;;;  
  ;; ;;   ;;;   ;; ;;     ;;    ;;    ;;     ;; ;;       ;;     ;; 
 ;;   ;;  ;;;;  ;; ;;     ;;    ;;    ;;     ;; ;;       ;;     ;; 
;;     ;; ;; ;; ;; ;;     ;;    ;;    ;;;;;;;;; ;;;;;;   ;;;;;;;;  
;;;;;;;;; ;;  ;;;; ;;     ;;    ;;    ;;     ;; ;;       ;;   ;;   
;;     ;; ;;   ;;; ;;     ;;    ;;    ;;     ;; ;;       ;;    ;;  
;;     ;; ;;    ;;  ;;;;;;;     ;;    ;;     ;; ;;;;;;;; ;;     ;; 

;;    ;; ;;;;;;;; ;;     ;; ;;;;;;;;   ;;;;;;;  ;;    ;; 
;;;   ;; ;;       ;;     ;; ;;     ;; ;;     ;; ;;;   ;; 
;;;;  ;; ;;       ;;     ;; ;;     ;; ;;     ;; ;;;;  ;; 
;; ;; ;; ;;;;;;   ;;     ;; ;;;;;;;;  ;;     ;; ;; ;; ;; 
;;  ;;;; ;;       ;;     ;; ;;   ;;   ;;     ;; ;;  ;;;; 
;;   ;;; ;;       ;;     ;; ;;    ;;  ;;     ;; ;;   ;;; 
;;    ;; ;;;;;;;;  ;;;;;;;  ;;     ;;  ;;;;;;;  ;;    ;; 

% show filters 
% compute filters from this neuron
clear K1 K2
for j = 1:size(data(1).S,2)
	S = data(1).S(:,j);
	X = data(1).X(:,j);
	R = data(1).R(:,j);

	temp1 = fitFilter2Data(S,X,'offset',200);
	temp2 = fitFilter2Data(S,R,'offset',200);
	K1(:,j) = temp1(100:end-100);
	K2(:,j) = temp2(100:end-100);
end
filtertime = 1e-3*(1:length(K1)) - .1;



% normalise the filters for display
nK1 = K1; nK2 = K2;
for i = 1:size(data(1).S,2)
	nK1(:,i) = K1(:,i)/norm(K1(:,i));
	nK2(:,i) = K2(:,i)/norm(K2(:,i));
end
nK2(:,isnan(sum(nK2))) = [];

filtertime = 1e-3*(1:length(nK1)) - .1;
errorShade(ax(5),filtertime,mean(nK1,2),sem(nK1'),'Color','k');
errorShade(ax(6),filtertime,mean(nK2,2),sem(nK2'),'Color','k');


% project using these filters

c = parula(6); 
c = c([1 2 4 5],:);
c = flipud(c);
ss = 10;

i = 1;
clear Xl Rl XL RL
for j = 1:size(data(i).S,2)
	S = data(i).S(:,j);
	X = data(i).X(:,j);
	R = data(i).R(:,j);

	Xp = convolve(time,S,K1(:,j),filtertime);
	K1(:,j) = K1(:,j)/(nanstd(Xp)/nanstd(S)); % normalise correctly 
	Xp = convolve(time,S,K1(:,j),filtertime);
	Xl(j) = plot(ax(7),Xp(1:ss:end),X(1:ss:end),'Color',c(j,:));
	XL{j} = ['r^2 = ' oval(rsquare(Xp,X))];

	if ~any(isnan(R))
		Rp = convolve(time,S,K2(:,j),filtertime);
		K2(:,j) = K2(:,j)/(nanstd(Rp)/nanstd(S)); % normalise correctly 
		Rp = convolve(time,S,K2(:,j),filtertime);
		Rl(j) = plot(ax(8),Rp(1:ss:end),R(1:ss:end),'Color',c(j,:));
		RL{j} = ['r^2 = ' oval(rsquare(Rp,R))];
	end

end
legend(Xl,XL,'Location','southeast')
legend(Rl(2:end),RL(2:end),'Location','southeast')



   ;;;    ;;;;;;;;   ;;;;;;;     ;;;    
  ;; ;;   ;;     ;; ;;     ;;   ;; ;;   
 ;;   ;;  ;;     ;;        ;;  ;;   ;;  
;;     ;; ;;;;;;;;   ;;;;;;;  ;;     ;; 
;;;;;;;;; ;;     ;;        ;; ;;;;;;;;; 
;;     ;; ;;     ;; ;;     ;; ;;     ;; 
;;     ;; ;;;;;;;;   ;;;;;;;  ;;     ;; 


cdata = consolidateData2(getPath(dataManager,'c2bce18a6b0a7e89e9c6832dcc27e39b'));
[cdata, data] =  assembleScaledNatStim(cdata);


% show filters 
% compute filters from this neuron
clear K1 K2
for j = 1:size(data(1).S,2)
	S = data(1).S(:,j);
	X = data(1).X(:,j);
	R = data(1).R(:,j);

	temp1 = fitFilter2Data(S,X,'offset',200);
	temp2 = fitFilter2Data(S,R,'offset',200);
	K1(:,j) = temp1(100:end-100);
	K2(:,j) = temp2(100:end-100);
end
filtertime = 1e-3*(1:length(K1)) - .1;

% normalise the filters for display
nK1 = K1; nK2 = K2;
for i = 1:size(data(1).S,2)
	nK1(:,i) = K1(:,i)/norm(K1(:,i));
	nK2(:,i) = K2(:,i)/norm(K2(:,i));
end
nK2(:,isnan(sum(nK2))) = [];

filtertime = 1e-3*(1:length(nK1)) - .1;
errorShade(ax(9),filtertime,mean(nK1,2),sem(nK1'),'Color','k');
errorShade(ax(10),filtertime,mean(nK2,2),sem(nK2'),'Color','k');


% project using these filters
c = parula(6); 
c = c([1 5],:);
c = flipud(c);
ss = 10;

i = 1;
clear Xl Rl XL RL
for j = 1:size(data(i).S,2)
	S = data(i).S(:,j);
	X = data(i).X(:,j);
	R = data(i).R(:,j);

	Xp = convolve(time,S,K1(:,j),filtertime);
	K1(:,j) = K1(:,j)/(nanstd(Xp)/nanstd(S)); % normalise correctly 
	Xp = convolve(time,S,K1(:,j),filtertime);
	Xl(j) = plot(ax(11),Xp(1:ss:end),X(1:ss:end),'Color',c(j,:));
	XL{j} = ['r^2 = ' oval(rsquare(Xp,X))];

	if ~any(isnan(R))
		Rp = convolve(time,S,K2(:,j),filtertime);
		K2(:,j) = K2(:,j)/(nanstd(Rp)/nanstd(S)); % normalise correctly 
		Rp = convolve(time,S,K2(:,j),filtertime);
		Rl(j) = plot(ax(12),Rp(1:ss:end),R(1:ss:end),'Color',c(j,:));
		RL{j} = ['r^2 = ' oval(rsquare(Rp,R))];
	end

end
legend(Xl,XL,'Location','southeast')
legend(Rl(1:end),RL(1:end),'Location','southeast')



;;          ;;;    ;;;;;;;;  ;;;;;;;; ;;        ;;;;;;          ;;;;;;;; ;;;;;;;;  ;;;;;;  
;;         ;; ;;   ;;     ;; ;;       ;;       ;;    ;;         ;;          ;;    ;;    ;; 
;;        ;;   ;;  ;;     ;; ;;       ;;       ;;               ;;          ;;    ;;       
;;       ;;     ;; ;;;;;;;;  ;;;;;;   ;;        ;;;;;;  ;;;;    ;;;;;;      ;;    ;;       
;;       ;;;;;;;;; ;;     ;; ;;       ;;             ;; ;;;;    ;;          ;;    ;;       
;;       ;;     ;; ;;     ;; ;;       ;;       ;;    ;;  ;;     ;;          ;;    ;;    ;; 
;;;;;;;; ;;     ;; ;;;;;;;;  ;;;;;;;; ;;;;;;;;  ;;;;;;  ;;      ;;;;;;;;    ;;     ;;;;;;  


% labels, prettification, etc
xlabel(ax(1),'Filter lag (s)')
ylabel(ax(1),'Filter (norm)')
title(ax(1),'ab2 LFP')
clear l L
l(1) = plot(ax(1),NaN,NaN,'k.','MarkerSize',24);
l(2) = plot(ax(1),NaN,NaN,'r.','MarkerSize',24);
L{1} = 'Nat. stimulus'; L{2} = 'Gaussian';
legend(l,L,'Location','southeast')

xlabel(ax(2),'Filter lag (s)')
ylabel(ax(2),'Filter (norm)')
title(ax(2),'ab2A firing rate')
clear l L
l(1) = plot(ax(2),NaN,NaN,'k.','MarkerSize',24);
l(2) = plot(ax(2),NaN,NaN,'r.','MarkerSize',24);
L{1} = 'Nat. stimulus'; L{2} = 'Gaussian';
legend(l,L,'Location','northeast')

xlabel(ax(3),['Stimulus projected with' char(10) 'Gaussian filter (V)'])
ylabel(ax(3),'ab2 LFP (mV)')
set(ax(3),'YDir','reverse','XDir','reverse')

xlabel(ax(4),['Stimulus projected with' char(10) 'Gaussian filter (V)'])
ylabel(ax(4),'Firing rate (Hz)')

xlabel(ax(5),'Filter lag (s)')
ylabel(ax(5),'Filter (norm)')
title(ax(5),'ab2 LFP')

xlabel(ax(6),'Filter lag (s)')
ylabel(ax(6),'Filter (norm)')
title(ax(6),'ab2A firing rate')

xlabel(ax(7),'Projected Stimulus (mV)')
ylabel(ax(7),'ab2 LFP (mV)')
set(ax(7),'XDir','reverse','YDir','reverse')

xlabel(ax(8),'Projected Stimulus (mV)')
ylabel(ax(8),'Firing rate (Hz)')


xlabel(ax(9),'Filter lag (s)')
ylabel(ax(9),'Filter (norm)')
title(ax(9),'ab3 LFP')

xlabel(ax(10),'Filter lag (s)')
ylabel(ax(10),'Filter (norm)')
title(ax(10),'ab3A firing rate')

xlabel(ax(11),'Projected Stimulus (mV)')
ylabel(ax(11),'ab3 LFP (mV)')
set(ax(11),'XDir','reverse','YDir','reverse')

xlabel(ax(12),'Projected Stimulus (mV)')
ylabel(ax(12),'ab3A firing rate (Hz)')


prettyFig();

if being_published
	snapnow
	delete(gcf)
end

clearvars -except being_published


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
end


figure('outerposition',[0 0 1403 801],'PaperUnits','points','PaperSize',[1403 801]); hold on
% first show the best fit model
x = logspace(-2,2,100);
y = 1./(1+(p.k_D./x).^p.n);
subplot(2,4,1); hold on
plot(x,y,'k')
set(gca,'XScale','log')
xlabel('Stimulus')
ylabel('a')


subplot(2,4,5); hold on
filtertime = 1:length(K); filtertime = filtertime - 50;
errorShade(filtertime,mean(K,2),std(K,[],2),'Color','k');
xlabel('Filter lag (ms)')
ylabel('Filter (norm)')
set(gca,'YDir','reverse')

% show one trace
example_trace = 3;
subplot(2,4,2:4); hold on
time = 1:length(data(2).S); time = time*1e-3;
clear l
l(1) = plot(time,data(2).X(:,example_trace),'k');
l(2) = plot(time,data(2).XP(:,example_trace),'r');
set(gca,'XLim',[0 70],'YDir','reverse')
xlabel('Time (s)')
ylabel('\Delta LFP (mV)')
r2 = rsquare(data(2).X(:,example_trace),data(2).XP(:,example_trace));
legend(l,{'ab2 LFP',['NL model, r^2 = ' oval(r2)]},'Location','northwest')


ax(6) = subplot(2,4,6); hold on
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

		plot(ax(6),t,S(a:z))
		plot(ax(7),t,X(a:z))
		plot(ax(8),t,XP(a:z))

	end
end

xlabel(ax(6),'Time since whiff (ms)')
ylabel(ax(6),'Stimulus (V)')

xlabel(ax(7),'Time since whiff (ms)')
ylabel(ax(7),'ab2 LFP (mV)')
set(ax(7),'YDir','reverse')

xlabel(ax(8),'Time since whiff (ms)')
ylabel(ax(8),'NL model (mV)')
set(ax(8),'YDir','reverse')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end


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

if being_published
	snapnow
	delete(gcf)
end

%% Version Info
%
pFooter;


