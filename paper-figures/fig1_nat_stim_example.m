% fig1_nat_stim_example

%%
% This script makes fig 1 for the eLife resubmission. The point of this figure is to show our naturalistic stimulus data, and the features of the ORN response we see in it. 

pHeader;


% make the figure and all subplots 
figure('outerposition',[0 0 1302 1001],'PaperUnits','points','PaperSize',[1000 1001]); hold on
clear ax
ax(1) = subplot(4,4,1:3); hold on
ax(2) = subplot(4,4,5:7); hold on
ax(3) = subplot(4,4,9:11); hold on

ax(4) = subplot(4,4,13); hold on
ax(5) = subplot(4,4,14); hold on

ax(6) = subplot(4,4,4); hold on
ax(7) = subplot(4,4,8); hold on
ax(8) = subplot(4,4,12); hold on


ax(9) = subplot(4,4,15); hold on
ax(10) = subplot(4,4,16); hold on


% get the data

cdata = consolidateData2(getPath(dataManager,'4608c42b12191b383c84fea52392ea97'));
[cdata, data] =  assembleScaledNatStim(cdata);
time = 1e-3*(1:length(data(1).S));

% plot the data
c = parula(6); 
c = c([1 3 5],:);
c = flipud(c);
for i = 1:3
	plot(ax(1),time,data(2).S(:,i),'Color',c(i,:))
	plot(ax(2),time,data(2).X(:,i),'Color',c(i,:))
	plot(ax(3),time,data(2).R(:,i),'Color',c(i,:))
end

% show whiff statistics 
i = 2;
for j = 1:size(data(i).S,2)
	S = data(i).S(:,j);
	X = data(i).X(:,j);
	R = data(i).R(:,j);
	ws = whiffStatistics(S,X,R,300,'MinPeakProminence',max(S/1e2),'debug',false);
	plot(ax(4),ws.stim_peaks,ws.peak_LFP,'.','MarkerSize',20,'Color',c(j,:))
	plot(ax(5),ws.stim_peaks,ws.peak_firing_rate,'.','MarkerSize',20,'Color',c(j,:))
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

		plot(ax(6),t,S(a:z))
		plot(ax(7),t,X(a:z))
		plot(ax(8),t,R(a:z))

	end
end

% show linear analysis output plot


ss = 50;

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



prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

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


%% Version Info
%
pFooter;

