
pHeader;

%% Analysis of naturalistic stimulus responses
% In this document I analyse all the data we have where neurons respond to naturalistic odorant stimuli. 


   ;;;    ;;;;;;;;   ;;;;;;;      ;;;;;;;     ;;;     ;;;;;;  
  ;; ;;   ;;     ;; ;;     ;;    ;;     ;;   ;; ;;   ;;    ;; 
 ;;   ;;  ;;     ;;        ;;           ;;  ;;   ;;  ;;       
;;     ;; ;;;;;;;;   ;;;;;;;      ;;;;;;;  ;;     ;; ;;       
;;;;;;;;; ;;     ;;        ;;    ;;        ;;;;;;;;; ;;       
;;     ;; ;;     ;; ;;     ;;    ;;        ;;     ;; ;;    ;; 
;;     ;; ;;;;;;;;   ;;;;;;;     ;;;;;;;;; ;;     ;;  ;;;;;;  


%% ab3A -- ethyl acetate
% The following figure shows the stimulus and responses from ab3A neurons to ethyl acetate odourant. 

% load the data
cdata = consolidateKontrollerData(getPath(dataManager,'f70e37a7db469b88c0fc79ff5e828e9d'));
v2struct(cdata); clear cdata



% remove baseline from the stimulus and the LFP
for i = 1:size(LFP,2)
	LFP(:,i) = LFP(:,i) - mean(LFP(1:5e3,i));
	PID(:,i) = PID(:,i) - min(PID(1:5e3,i));
end

time = 1e-3*(1:length(PID));

figure('outerposition',[0 0 1000 901],'PaperUnits','points','PaperSize',[1000 901]); hold on
subplot(3,1,1); hold on
c = lines(max(orn));
for i = max(orn):-1:1
	S = nanmean(PID(:,orn==i),2);
	plot(time,S,'Color',c(i,:))
end
ylabel('Stimulus (V)')
set(gca,'XLim',[0 70])

subplot(3,1,2); hold on
for i = max(orn):-1:1
	X = nanmean(LFP(:,orn==i),2);
	plot(time,X,'Color',c(i,:))
end
ylabel('\Delta LFP (mV)')
set(gca,'XLim',[0 70])

subplot(3,1,3); hold on
for i = max(orn):-1:1
	X = fA(:,orn == i);
	X(:,sum(X)==0) = [];
	if ~isempty(X)
		plot(time,nanmean(X,2),'Color',c(i,:))
	end
end
ylabel('Firing rate (Hz)')
xlabel('Time (s)')
set(gca,'XLim',[0 70])

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end


%%
% How does the neuron respond to every whiff of the stimulus? How does the neuron handle the very broad distribution of the stimulus, given that its response range is limited? In the following figure, I plot, for each neuron, for each whiff, the maximum stimulus and maximum firing rate response. 

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
for i = 1:max(orn)
	X = fA(:,orn == i);
	X(:,sum(X)==0) = [];
	R = nanmean(X,2);
	S = nanmean(PID(:,orn==i),2);
	ws = whiffStatistics(S,X,R,300);
	plot(ws.stim_peaks,ws.peak_firing_rate,'.','MarkerSize',20,'Color',c(i,:))
end
set(gca,'XScale','log','YScale','log','XLim',[.1 10])
xlabel('S_{max} (V)')
ylabel('R_{max} (Hz)')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end



%%
% Now, I find all whiffs in the stimulus, and plot them, coloured by peak. I also plot the corresponding LFP and firing rate responses. 

figure('outerposition',[0 0 1523 901],'PaperUnits','points','PaperSize',[1523 901]); hold on
norns = 6;
for i = 1:norns
	
	X = fA(:,orn == i);
	X(:,sum(X)==0) = [];
	R = nanmean(X,2);
	X = nanmean(LFP(:,orn==i),2);
	S = nanmean(PID(:,orn==i),2);
	ws = whiffStatistics(S,X,R,300);
	[ws.stim_peaks,idx] = sort(ws.stim_peaks);
	ws.stim_peak_loc = ws.stim_peak_loc(idx);
	c = parula(length(idx));

	a = ws.stim_peak_loc - 200;
	z = ws.stim_peak_loc + 200;

	subplot(3,norns,i); hold on
	for j = 1:length(idx)
		plot(S(a(j):z(j)),'Color',c(j,:))
	end
	ylabel('Stimulus (V)')

	subplot(3,norns,i+norns); hold on
	for j = 1:length(idx)
		plot(X(a(j):z(j)),'Color',c(j,:))
	end
	ylabel('\Delta LFP (mV)')

	subplot(3,norns,i+2*norns); hold on
	for j = 1:length(idx)
		plot(R(a(j):z(j)),'Color',c(j,:))
	end
	ylabel('Firing rate (Hz)')

end

prettyFig();

if being_published
	snapnow
	delete(gcf)
end



   ;;;    ;;;;;;;;   ;;;;;;;  
  ;; ;;   ;;     ;; ;;     ;; 
 ;;   ;;  ;;     ;;        ;; 
;;     ;; ;;;;;;;;   ;;;;;;;  
;;;;;;;;; ;;     ;;        ;; 
;;     ;; ;;     ;; ;;     ;; 
;;     ;; ;;;;;;;;   ;;;;;;;  


 ;;;;;;;          ;;;;;;;;  ;;     ;; ;;;;;;;;    ;;;    ;;    ;;  ;;;;;;;  ;;    ;; ;;;;;;;; 
;;     ;;         ;;     ;; ;;     ;;    ;;      ;; ;;   ;;;   ;; ;;     ;; ;;;   ;; ;;       
       ;;         ;;     ;; ;;     ;;    ;;     ;;   ;;  ;;;;  ;; ;;     ;; ;;;;  ;; ;;       
 ;;;;;;;  ;;;;;;; ;;;;;;;;  ;;     ;;    ;;    ;;     ;; ;; ;; ;; ;;     ;; ;; ;; ;; ;;;;;;   
;;                ;;     ;; ;;     ;;    ;;    ;;;;;;;;; ;;  ;;;; ;;     ;; ;;  ;;;; ;;       
;;                ;;     ;; ;;     ;;    ;;    ;;     ;; ;;   ;;; ;;     ;; ;;   ;;; ;;       
;;;;;;;;;         ;;;;;;;;   ;;;;;;;     ;;    ;;     ;; ;;    ;;  ;;;;;;;  ;;    ;; ;;;;;;;; 


%% ab3A -- 2-butanone
% The following figure shows the stimulus and responses from ab3A neurons to 2-butanone odourant. 

% load the data
cdata = consolidateKontrollerData(getPath(dataManager,'b8d40a4b987ccd1926bbd6d4578bbd99'));
v2struct(cdata); clear cdata


% remove baseline from the stimulus and the LFP
for i = 1:size(LFP,2)
	LFP(:,i) = LFP(:,i) - mean(LFP(1:5e3,i));
	PID(:,i) = PID(:,i) - min(PID(1:5e3,i));
end

time = 1e-3*(1:length(PID));

figure('outerposition',[0 0 1000 901],'PaperUnits','points','PaperSize',[1000 901]); hold on
subplot(3,1,1); hold on
c = lines(max(orn));
for i = max(orn):-1:1
	S = nanmean(PID(:,orn==i),2);
	plot(time,S,'Color',c(i,:))
end
ylabel('Stimulus (V)')
set(gca,'XLim',[0 70])

subplot(3,1,2); hold on
for i = max(orn):-1:1
	X = nanmean(LFP(:,orn==i),2);
	plot(time,X,'Color',c(i,:))
end
ylabel('\Delta LFP (mV)')
set(gca,'XLim',[0 70])

subplot(3,1,3); hold on
for i = max(orn):-1:1
	X = fA(:,orn == i);
	X(:,sum(X)==0) = [];
	if ~isempty(X)
		plot(time,nanmean(X,2),'Color',c(i,:))
	end
end
ylabel('Firing rate (Hz)')
xlabel('Time (s)')
set(gca,'XLim',[0 70])

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%%
% How does the neuron respond to every whiff of the stimulus? How does the neuron handle the very broad distribution of the stimulus, given that its response range is limited? In the following figure, I plot, for each neuron, for each whiff, the maximum stimulus and maximum firing rate response. 

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
for i = 1:2
	X = fA(:,orn == i);
	X(:,sum(X)==0) = [];
	R = nanmean(X,2);
	S = nanmean(PID(:,orn==i),2);
	ws = whiffStatistics(S,X,R,300);
	plot(ws.stim_peaks,ws.peak_firing_rate,'.','MarkerSize',20,'Color',c(i,:))
end
set(gca,'XScale','log','YScale','log')
xlabel('S_{max} (V)')
ylabel('R_{max} (Hz)')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% Now, I find all whiffs in the stimulus, and plot them, coloured by peak. I also plot the corresponding LFP and firing rate responses. 

figure('outerposition',[0 0 1200 901],'PaperUnits','points','PaperSize',[1200 901]); hold on
for i = 1:2
	
	X = fA(:,orn == i);
	X(:,sum(X)==0) = [];
	R = nanmean(X,2);
	X = nanmean(LFP(:,orn==i),2);
	S = nanmean(PID(:,orn==i),2);
	ws = whiffStatistics(S,X,R,300);
	[ws.stim_peaks,idx] = sort(ws.stim_peaks);
	ws.stim_peak_loc = ws.stim_peak_loc(idx);
	c = parula(length(idx));

	a = ws.stim_peak_loc - 200;
	z = ws.stim_peak_loc + 200;

	subplot(2,3,(i-1)*3+1); hold on
	for j = 1:length(idx)
		plot(S(a(j):z(j)),'Color',c(j,:))
	end
	ylabel('Stimulus (V)')

	subplot(2,3,(i-1)*3+2); hold on
	for j = 1:length(idx)
		plot(X(a(j):z(j)),'Color',c(j,:))
	end
	ylabel('\delta LFP (mV)')

	subplot(2,3,(i-1)*3+3); hold on
	for j = 1:length(idx)
		plot(R(a(j):z(j)),'Color',c(j,:))
	end
	ylabel('Firing rate (Hz)')

end

prettyFig();

if being_published
	snapnow
	delete(gcf)
end


%%
% Now, I plot the firing rate responses as a function of the linear projection for each neuron 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
for i = 1:2
	subplot(1,2,i); hold on
	X = fA(:,orn == i);
	X(:,sum(X)==0) = [];
	R = nanmean(X,2);
	S = nanmean(PID(:,orn==i),2);
	K = fitFilter2Data(S,R,'filter_length',1e3,'offset',200);
	K = K(100:end-100); filtertime = 1e-3*(1:length(K)) - .1;
	fp = convolve(time,S,K,filtertime);
	fp = fp*nanstd(S)/nanstd(fp);
	plot(fp,R,'k.')
	xlabel('Projected stimulus (V)')
	ylabel('ab3A firing rate (Hz)')
end

prettyFig();

if being_published
	snapnow
	delete(gcf)
end


   ;;;    ;;;;;;;;   ;;;;;;;     ;;;        ;;;;;;;     ;;;     ;;;;;;  
  ;; ;;   ;;     ;; ;;     ;;   ;; ;;      ;;     ;;   ;; ;;   ;;    ;; 
 ;;   ;;  ;;     ;;        ;;  ;;   ;;            ;;  ;;   ;;  ;;       
;;     ;; ;;;;;;;;   ;;;;;;;  ;;     ;;     ;;;;;;;  ;;     ;; ;;       
;;;;;;;;; ;;     ;;        ;; ;;;;;;;;;    ;;        ;;;;;;;;; ;;       
;;     ;; ;;     ;; ;;     ;; ;;     ;;    ;;        ;;     ;; ;;    ;; 
;;     ;; ;;;;;;;;   ;;;;;;;  ;;     ;;    ;;;;;;;;; ;;     ;;  ;;;;;;  

;;;;;;;;  ;;;;;;;; ;;    ;;  ;;;;;;  ;;;;;;;; 
;;     ;; ;;       ;;;   ;; ;;    ;; ;;       
;;     ;; ;;       ;;;;  ;; ;;       ;;       
;;     ;; ;;;;;;   ;; ;; ;;  ;;;;;;  ;;;;;;   
;;     ;; ;;       ;;  ;;;;       ;; ;;       
;;     ;; ;;       ;;   ;;; ;;    ;; ;;       
;;;;;;;;  ;;;;;;;; ;;    ;;  ;;;;;;  ;;;;;;;; 

%% Dense naturalistic stimulus
% Now, I plot the per-whiff responses of ab3A neurons to ethyl acetate stimulus in the dense naturalistic stimulus case. 

clear od
load(getPath(dataManager,'aeb361c027b71938021c12a6a12a85cd'),'-mat')

figure('outerposition',[0 0 600 601],'PaperUnits','points','PaperSize',[600 601]); hold on
c = lines(length(od));
for i = 1:length(od)
	S = nanmean(od(i).stimulus,2); S = S - mean(S(1:5e3));
	R = od(i).firing_rate; 
	R(:,max(R) == 0 | isnan(max(R))) = [];
	R = nanmean(R,2);
	ws = whiffStatistics(S,R,R,300,'MinPeakProminence',max(S/1e2),'debug',false);
	plot(ws.stim_peaks,ws.peak_firing_rate,'.','MarkerSize',20,'Color',c(i,:))
	set(gca,'YLim',[0 150],'XScale','log')
	xlabel('Stim_{peak} (V)')
	ylabel('ab3A resp_{peak} (Hz)')
end

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

 ;;;;;;   ;;;;;;     ;;;    ;;       ;;;;;;;; ;;;;;;;;  
;;    ;; ;;    ;;   ;; ;;   ;;       ;;       ;;     ;; 
;;       ;;        ;;   ;;  ;;       ;;       ;;     ;; 
 ;;;;;;  ;;       ;;     ;; ;;       ;;;;;;   ;;     ;; 
      ;; ;;       ;;;;;;;;; ;;       ;;       ;;     ;; 
;;    ;; ;;    ;; ;;     ;; ;;       ;;       ;;     ;; 
 ;;;;;;   ;;;;;;  ;;     ;; ;;;;;;;; ;;;;;;;; ;;;;;;;;  

;;    ;;    ;;;    ;;;;;;;;     ;;;;;;  ;;;;;;;; ;;;; ;;     ;; 
;;;   ;;   ;; ;;      ;;       ;;    ;;    ;;     ;;  ;;;   ;;; 
;;;;  ;;  ;;   ;;     ;;       ;;          ;;     ;;  ;;;; ;;;; 
;; ;; ;; ;;     ;;    ;;        ;;;;;;     ;;     ;;  ;; ;;; ;; 
;;  ;;;; ;;;;;;;;;    ;;             ;;    ;;     ;;  ;;     ;; 
;;   ;;; ;;     ;;    ;;       ;;    ;;    ;;     ;;  ;;     ;; 
;;    ;; ;;     ;;    ;;        ;;;;;;     ;;    ;;;; ;;     ;; 

%% Scaled Naturalisitc Stimulus: ab2A and 2-butanone 
% In this section, I plot responses of ORNs to the naturalistic stimulus and scaled versions of the naturalistic stimulus. Each data set is from the same neuron. 


   ;;;    ;;;;;;;;   ;;;;;;;           ;;;;;;;  ;;;;;;;;  ;;     ;; ;;;;;;;; 
  ;; ;;   ;;     ;; ;;     ;;         ;;     ;; ;;     ;; ;;     ;;    ;;    
 ;;   ;;  ;;     ;;        ;;                ;; ;;     ;; ;;     ;;    ;;    
;;     ;; ;;;;;;;;   ;;;;;;;  ;;;;;;;  ;;;;;;;  ;;;;;;;;  ;;     ;;    ;;    
;;;;;;;;; ;;     ;; ;;                ;;        ;;     ;; ;;     ;;    ;;    
;;     ;; ;;     ;; ;;                ;;        ;;     ;; ;;     ;;    ;;    
;;     ;; ;;;;;;;;  ;;;;;;;;;         ;;;;;;;;; ;;;;;;;;   ;;;;;;;     ;;    


% get all data 
cdata = consolidateData2(getPath(dataManager,'4608c42b12191b383c84fea52392ea97'));
[cdata, data] =  assembleScaledNatStim(cdata);

t = 'ab2A -- 2-butanone';
plotScaledNatStimData(data,being_published,t);

%%
% Now, I plot the peak LFP and firing rate response in each whiff as a function of whiff intensity, across all the stimulus presentations. 

plotScaledNatStimWhiffStats(data);
suptitle(t)
prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% From the previous figure, it looks as though there are multiple responses for whiffs of the same magnitude. Is this true? To look into this more closely, I collect whiffs with the same amplitude, and plot the corresponding LFP and firing rate responses. I see that there is variation in the LFP and firing rate responses, and this variation seems to be explained by responses to preceding stimuli. 

s_range = [.8 1.2];

figure('outerposition',[0 0 1501 502],'PaperUnits','points','PaperSize',[1501 502]); hold on
for i = 1:3
	ax(i) = subplot(1,3,i); hold on
end

show_these = [           2       27764
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

		plot(ax(1),S(a:z))
		plot(ax(2),X(a:z))
		plot(ax(3),R(a:z))

	end
end
ylabel(ax(1),'Stimulus (V)')
ylabel(ax(2),'LFP (mV)')
ylabel(ax(3),'Firing rate (Hz)')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% Now, I extract filters for each case and plot the filters and linear projectins on a per-neuron basis, for the LFP.

plotScaledNatStimFilters(data,'LFP')
suptitle(t)
prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% Now, I extract filters for each case and plot the filters and linear projectins on a per-neuron basis, for the LFP.

plotScaledNatStimFilters(data,'firing')
suptitle(t)
prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% Can we account for the firing rate responses of this ORN using a NLN model? 

this_orn = 2;
clear fd
for i = 1:size(data(this_orn).X,2)
	S = data(this_orn).S(:,i); S = S - min(S);
	R = data(this_orn).R(:,i);
	fd(i).stimulus = [S R];
	fd(i).response = R;
	%fd(i).response(fd(i).response<10) = NaN;
end

clear p
p.  k_D = 0.0933;
p.    n = 1;

% generate responses using this model 
for j = 1:size(data(2).X,2)
	[data(2).P(:,j),K(:,j)] = NLNmodel([data(2).S(:,j) - min(data(2).S(:,j)) data(2).R(:,j)] ,p);
end
nK = K;
for i = 1:3
	temp = fitFilter2Data(data(2).S(:,i),data(2).R(:,i),'filter_length',700,'offset',100);
	nK(:,i) = temp(50:end-50);
end

figure('outerposition',[0 0 1601 901],'PaperUnits','points','PaperSize',[1601 901]); hold on
subplot(2,4,1:3); hold on
time = 1e-3*(1:length(data(2).R));
plot(time,data(2).R(:,2),'k')
plot(time,data(2).P(:,2),'r')
legend({'ab2A','NLN model'})
xlabel('Time (s)')
ylabel('Firing rate (Hz)')

% compare nonlinearity to what we saw from the whiffs
subplot(2,4,5); hold on
for i = 2
	for j = 1:3
		S = data(i).S(:,j);
		R = data(i).R(:,j);
		ws = whiffStatistics(S,R,R,300,'MinPeakProminence',max(S/1e2),'debug',false);
		plot(ws.stim_peaks,ws.peak_firing_rate,'.','MarkerSize',20,'Color','k')
	end
end
set(gca,'XScale','log')
% plot the input nonlinearity
x = logspace(-4,2,100);
a = 1./(1+(p.k_D./x).^p.n);
[ax,~,l] = plotyy(x,NaN*x,x,a);
ax(2).XScale = 'log';
ax(1).YLim = [0 250];
set(ax,'XLim',[1e-3 1e2],'XTick',[1e-3 1e-2 1e-1 1e0 1e1 1e2])
xlabel('Whiff intensity (V)')
ylabel(ax(1),'ORN response (Hz)')
ylabel(ax(2),'Input nonlinearity')

% compare the filter to the normal filter
subplot(2,4,6); hold on
% extract normal filters

for i = 1:3
	nK(:,i) = nK(:,i)/norm(nK(:,i));
	K(:,i) = K(:,i)/norm(K(:,i));
end
filtertime = 1e-3*(1:length(K)) - .05;
errorShade(filtertime,mean(nK,2),sem(K'),'Color','k');
errorShade(filtertime,mean(K,2),sem(K'),'Color','r');
clear l
l(1) = plot(NaN,NaN,'k');
l(2) = plot(NaN,NaN,'r');
legend(l,{'Filter','with input NL'})
xlabel('Filter lag (s)')


% compare per-whiff responses across everything
subplot(2,4,4); hold on
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
plot(x,y,'k.','MarkerSize',20)
legend(['r^2 = ' oval(rsquare(x,y))],'Location','southeast')
xlabel('NLN model response (Hz)')
ylabel('ab2A response (Hz)')
title('Whiff-specific responses')
set(gca,'XLim',[0 300],'YLim',[0 300])

% now make a plot of the NLN responses vs. linear projections thereof, and show that we get the same sort of structure. 

subplot(2,4,7); hold on
for i = 1:3
	S = data(2).S(:,i);
	P = data(2).P(:,i);
	temp = fitFilter2Data(S,P,'filter_length',700,'offset',100);
	temp = temp(50:end-50);
	fp = convolve(time,S,temp,filtertime);
	plot(fp,P,'.')
end
xlabel('Proj. Stim. (V)')
ylabel('NLN firing rate (Hz)')


prettyFig();

if being_published
	snapnow
	delete(gcf)
end

s_range = [.8 1.2];

figure('outerposition',[0 0 1501 502],'PaperUnits','points','PaperSize',[1501 502]); hold on
for i = 1:3
	ax(i) = subplot(1,3,i); hold on
end

show_these = [           2       27764
           2       28790
           2       59776
           3       58049];

for i = 2
	c = lines(length(show_these));
	for j = 1:length(show_these)
		this_stim = show_these(j,1);
		this_loc = show_these(j,2);

		S = data(i).S(:,this_stim);
		X = data(i).X(:,this_stim);
		R = data(i).R(:,this_stim);
		P = data(i).P(:,this_stim);

		a = this_loc - 300;
		z = this_loc+300;

		plot(ax(1),S(a:z))
		plot(ax(2),X(a:z))
		plot(ax(3),R(a:z),'Color',c(j,:))
		plot(ax(3),P(a:z),':','Color',c(j,:))

	end
end
ylabel(ax(1),'Stimulus (V)')
ylabel(ax(2),'LFP (mV)')
ylabel(ax(3),'Firing rate (Hz)')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end



   ;;;    ;;;;;;;;   ;;;;;;;     ;;;             ;;;;;;;  ;;;;;;;;  ;;     ;; ;;;;;;;; 
  ;; ;;   ;;     ;; ;;     ;;   ;; ;;           ;;     ;; ;;     ;; ;;     ;;    ;;    
 ;;   ;;  ;;     ;;        ;;  ;;   ;;                 ;; ;;     ;; ;;     ;;    ;;    
;;     ;; ;;;;;;;;   ;;;;;;;  ;;     ;; ;;;;;;;  ;;;;;;;  ;;;;;;;;  ;;     ;;    ;;    
;;;;;;;;; ;;     ;;        ;; ;;;;;;;;;         ;;        ;;     ;; ;;     ;;    ;;    
;;     ;; ;;     ;; ;;     ;; ;;     ;;         ;;        ;;     ;; ;;     ;;    ;;    
;;     ;; ;;;;;;;;   ;;;;;;;  ;;     ;;         ;;;;;;;;; ;;;;;;;;   ;;;;;;;     ;;    

%% ab3A and 2-butanone.
% Now I make the same plots, but for a ab3A neuron responding to 2-butanone. 

% get all data 
cdata = consolidateData2(getPath(dataManager,'c2bce18a6b0a7e89e9c6832dcc27e39b'));
[cdata, data] =  assembleScaledNatStim(cdata);

t = 'ab3A -- 2-butanone';
plotScaledNatStimData(data,being_published,t);

%%
% Now, for each neuron, I plot the compare the stimulus, LFP and firing rate for the various doses for each neuron. Since the rough temporal structure of the stimulus is the same, this lets us compare the dynamics and scale of the responses directly. In the following figure, each column is data from a single neuron. I plot the stimulus vs. its lowest variant, and the LFP vs. the LFP for the lowest stimulus variant, and so on. 

[fold_changes] = plotScaledNatStimFoldChange1(data,being_published);

%%
% Now, I compute the fold change in the stimulus, LFP and firing rate as the slope of the lines in these plots, and look at how scaling the stimulus up changes the response in the LFP and the firing rate. 


plotScaledNatStimFoldChange2(fold_changes,being_published);


%%
% Now, I plot the peak LFP and firing rate response in each whiff as a function of whiff intensity, across all the stimulus presentations. 

plotScaledNatStimWhiffStats(data);
suptitle(t)
prettyFig();

if being_published
	snapnow
	delete(gcf)
end


%%
% Now, I extract filters for each case and plot the filters and linear projectins on a per-neuron basis, for the LFP.

plotScaledNatStimFilters(data,'LFP')
suptitle(t)
prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% Now, I extract filters for each case and plot the filters and linear projectins on a per-neuron basis, for the LFP.

plotScaledNatStimFilters(data,'firing')
suptitle(t)
prettyFig();

if being_published
	snapnow
	delete(gcf)
end

   ;;;    ;;;;;;;;   ;;;;;;;     ;;;        ;;;;;;;     ;;;     ;;;;;;  
  ;; ;;   ;;     ;; ;;     ;;   ;; ;;      ;;     ;;   ;; ;;   ;;    ;; 
 ;;   ;;  ;;     ;;        ;;  ;;   ;;            ;;  ;;   ;;  ;;       
;;     ;; ;;;;;;;;   ;;;;;;;  ;;     ;;     ;;;;;;;  ;;     ;; ;;       
;;;;;;;;; ;;     ;;        ;; ;;;;;;;;;    ;;        ;;;;;;;;; ;;       
;;     ;; ;;     ;; ;;     ;; ;;     ;;    ;;        ;;     ;; ;;    ;; 
;;     ;; ;;;;;;;;   ;;;;;;;  ;;     ;;    ;;;;;;;;; ;;     ;;  ;;;;;;  


%% ab2A and ethyl acetate.
% Now I make the same plots, but for a ab2A neuron responding to ethyl acetate. 

% get all data 
cdata = consolidateData2(getPath(dataManager,'9f359fef5acd000e1a24902d33b460ee'));
[cdata, data] =  assembleScaledNatStim(cdata);

t = 'ab2A -- ethyl acetate';
plotScaledNatStimData(data,being_published,t);

%%
% Now, for each neuron, I plot the compare the stimulus, LFP and firing rate for the various doses for each neuron. Since the rough temporal structure of the stimulus is the same, this lets us compare the dynamics and scale of the responses directly. In the following figure, each column is data from a single neuron. I plot the stimulus vs. its lowest variant, and the LFP vs. the LFP for the lowest stimulus variant, and so on. 

[fold_changes] = plotScaledNatStimFoldChange1(data,being_published);

%%
% Now, I compute the fold change in the stimulus, LFP and firing rate as the slope of the lines in these plots, and look at how scaling the stimulus up changes the response in the LFP and the firing rate. 


plotScaledNatStimFoldChange2(fold_changes,being_published);



%%
% Now, I plot the peak LFP and firing rate response in each whiff as a function of whiff intensity, across all the stimulus presentations. 

plotScaledNatStimWhiffStats(data);
suptitle(t)
prettyFig();

if being_published
	snapnow
	delete(gcf)
end


%%
% Now, I extract filters for each case and plot the filters and linear projectins on a per-neuron basis, for the LFP.

plotScaledNatStimFilters(data,'LFP')
suptitle(t)
prettyFig();

if being_published
	snapnow
	delete(gcf)
end



%% Version Info
%
pFooter;

