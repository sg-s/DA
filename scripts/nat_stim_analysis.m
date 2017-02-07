
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

%
% Now, I find all whiffs in the stimulus, and plot them, coloured by peak. I also plot the corresponding LFP and firing rate responses. 

% figure('outerposition',[0 0 1200 901],'PaperUnits','points','PaperSize',[1200 901]); hold on
% for i = 1:2
	
% 	X = fA(:,orn == i);
% 	X(:,sum(X)==0) = [];
% 	R = nanmean(X,2);
% 	X = nanmean(LFP(:,orn==i),2);
% 	S = nanmean(PID(:,orn==i),2);
% 	ws = whiffStatistics(S,X,R,300);
% 	[ws.stim_peaks,idx] = sort(ws.stim_peaks);
% 	ws.stim_peak_loc = ws.stim_peak_loc(idx);
% 	c = parula(length(idx));

% 	a = ws.stim_peak_loc - 200;
% 	z = ws.stim_peak_loc + 200;

% 	subplot(2,3,(i-1)*3+1); hold on
% 	for j = 1:length(idx)
% 		plot(S(a(j):z(j)),'Color',c(j,:))
% 	end
% 	ylabel('Stimulus (V)')

% 	subplot(2,3,(i-1)*3+2); hold on
% 	for j = 1:length(idx)
% 		plot(X(a(j):z(j)),'Color',c(j,:))
% 	end
% 	ylabel('\delta LFP (mV)')

% 	subplot(2,3,(i-1)*3+3); hold on
% 	for j = 1:length(idx)
% 		plot(R(a(j):z(j)),'Color',c(j,:))
% 	end
% 	ylabel('Firing rate (Hz)')

% end

% prettyFig();

% if being_published
% 	snapnow
% 	delete(gcf)
% end


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

figure('outerposition',[0 0 600 601],'PaperUnits','points','PaperSize',[1200 601]); hold on
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
% Can I account for the LFP responses of the ORN using a NL model? In the following figure, I fit a NL model to the LFP responses, and compare model predictions to data. 

this_orn = 2;
clear fd
for i = 1:size(data(this_orn).X,2)
	S = data(this_orn).S(:,i); S = S - min(S);
	R = data(this_orn).X(:,i);
	%fd(i).stimulus = [S R];
	fd(i).stimulus = S;
	fd(i).response = -R;
	%fd(i).response(fd(i).response<10) = NaN;
end

clear p
p.k_D = 0.1996;
p.n = 1.5000;

% generate responses using this model 
for j = 1:size(data(2).X,2)
	[data(2).XP(:,j)] = NLModel([data(2).S(:,j) - min(data(2).S(:,j)) data(2).X(:,j)] ,p);
end


figure('outerposition',[0 0 1501 502],'PaperUnits','points','PaperSize',[1501 502]); hold on
for i = 1:3
	ax(i) = subplot(1,3,i); hold on
end

clear l L
for i = 1:3
	l(i) = plot(ax(1),data(2).XP(:,i),data(2).X(:,i));
	L{i} = ['r^2 = ' oval(rsquare(data(2).XP(:,i),data(2).X(:,i)))];
end
legend(l,L,'Location','northwest')
xlabel(ax(1),'NL model')
ylabel(ax(1),'\Delta LVP (mV)')

show_these = [2       27764
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
		P = data(i).XP(:,this_stim);

		a = this_loc - 300;
		z = this_loc+300;

		plot(ax(2),S(a:z))
		plot(ax(3),X(a:z),'Color',c(j,:))
		plot(ax(3),P(a:z),':','Color',c(j,:))

	end
end
ylabel(ax(2),'Stimulus (V)')
ylabel(ax(3),'LFP (mV)')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% Even though the non-adapting NL model seems to do a good job, it doesn't seem to be able to capture the context-dependent modulation of LFP responses. What if I fit an adapting NL model to this? Can the addition of an adapting $k_D$ account for the context-dependent response modulation we see in the LFP? In the following figure, I fit an adapting NL model to the LFP data, and examine the same whiffs to see if it has captured the observed variation in the data. 


clear p
p.  k0 = 0.2720;
p.   B = 10;
p.   n = 2.3750;
p.   A = 0.6992;
p.tauA = 13.8633;
p.tauB = 57.6250;
p.   C = 2.6069;
p.   D = -0.7344;



% generate responses using this model 
clear k_D
k_D = NaN*data(2).S;
for j = 1:size(data(2).X,2)
	[data(2).XP(:,j),~,~,k_D(:,j)] = aNL4(fd(j).stimulus ,p);
end
data(2).XP = -data(2).XP;


figure('outerposition',[0 0 1501 901],'PaperUnits','points','PaperSize',[1501 901]); hold on
ax(1) = subplot(2,3,1:3); hold on
for i = 1:3
	ax(i+1) = subplot(2,3,3+i); hold on
end
for i = 1:3
	plot(ax(1),1e-3*(1:length(k_D)),k_D(:,i))
end
xlabel(ax(1),'Time (s)')
ylabel(ax(1),'k_D')

clear l L
for i = 1:3
	l(i) = plot(ax(2),data(2).XP(:,i),data(2).X(:,i));
	L{i} = ['r^2 = ' oval(rsquare(data(2).XP(:,i),data(2).X(:,i)))];
end
legend(l,L,'Location','northwest')
xlabel(ax(2),'NL model')
ylabel(ax(2),'\Delta LFP (mV)')

show_these = [2       27764
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
		P = data(i).XP(:,this_stim);

		a = this_loc - 300;
		z = this_loc+300;

		plot(ax(3),S(a:z))
		plot(ax(4),X(a:z),'Color',c(j,:))
		plot(ax(4),P(a:z),':','Color',c(j,:))

	end
end
ylabel(ax(3),'Stimulus (V)')
ylabel(ax(4),'LFP (mV)')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% Now, I extract filters for each case and plot the filters and linear projectins on a per-neuron basis, for the firing rate.

plotScaledNatStimFilters(data,'firing')
suptitle(t)
prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% Can we account for the firing rate responses of this ORN using a NLN model? In the following figure, I fit a non-adapting NLN model to the firing rate data, and compare the model predictions (red) to the data (black). Then, I compare the responses of the model to the responses of the neuron for every single whiff, and find a good agreement. In the bottom row, I show the parameters of the NL model. Note that that the best-fit nonlinearity roughly follows the nonlinearity apparent in the whiff-specific responses. Note also that the filter computed in the presence of the nonlinearity is much sharper, and has a pronounced negative lobe. Finally, I plot the NLN firing rate as a function of the projected stimulus, to show that the NLN model has similar deviations from linearity that we see in the data. 


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
clear K
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

%%
% Can the NLN model account for the context-dependent variation we observe? I make the same pot as before where I collect whiffs of similar amplitude, and look at firing rate responses. I also plot the NLN model responses in dotted lines. It looks like the NLN model is partly able to account for this, presumable because of the negative lobe in the filter. 

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
% In this section I look at the responses of ab3A to 2-butanone. The following figure shows the stimulus, LFP and firing rate from a single ORN. 

% get all data 
cdata = consolidateData2(getPath(dataManager,'c2bce18a6b0a7e89e9c6832dcc27e39b'));
[cdata, data] =  assembleScaledNatStim(cdata);

t = 'ab3A -- 2-butanone';
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
% In this section, I plot the responses of ab2A to ethyl acetate. The following figure shows the stimulus, LFP and firing rate (where available -- the neuron pinches easily, so it's not always there). 

% get all data 
cdata = consolidateData2(getPath(dataManager,'9f359fef5acd000e1a24902d33b460ee'));
[cdata, data] =  assembleScaledNatStim(cdata);

t = 'ab2A -- ethyl acetate';
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


