pHeader;


%% One model to fit all the data 
% In this document I attempt to find one model that fits all the data we have so far. First, I'll fit some simplified version of the model to each data set, and study how well it constrains the parameters of this model. Then, finally, I'll attempt to fit one model to all the data, and compare that model's prediction to the prediction from models that are fit to one dataset alone. 

% first get all the data

 ;;;;;;   ;;;;;;;; ;;;;;;;;    ;;;;;;;;     ;;;    ;;;;;;;;    ;;;    
;;    ;;  ;;          ;;       ;;     ;;   ;; ;;      ;;      ;; ;;   
;;        ;;          ;;       ;;     ;;  ;;   ;;     ;;     ;;   ;;  
;;   ;;;; ;;;;;;      ;;       ;;     ;; ;;     ;;    ;;    ;;     ;; 
;;    ;;  ;;          ;;       ;;     ;; ;;;;;;;;;    ;;    ;;;;;;;;; 
;;    ;;  ;;          ;;       ;;     ;; ;;     ;;    ;;    ;;     ;; 
 ;;;;;;   ;;;;;;;;    ;;       ;;;;;;;;  ;;     ;;    ;;    ;;     ;; 


% get the filter from the Gaussian stimuli 
clear MSGdata
MSGdata = consolidateData2(getPath(dataManager,'93ba5d68174e3df9f462a1fc48c581da'));
MSGdata = cleanMSGdata(MSGdata);

% make sure stimulus is always positive
for i = 1:size(MSGdata.PID,2)
	MSGdata.PID(:,i) = MSGdata.PID(:,i) - min(MSGdata.PID(:,i));
end

% get the sparse naturalistic stimuli with LFP
% get this data
clearvars -except being_published MSGdata
load('/local-data/DA-paper/data-for-paper/nat-stim/ab3A_nat_stim.ORNData','-mat')
od(1) = [];
SNSdata.LFP = [];
SNSdata.fA = [];
SNSdata.PID = [];
SNSdata.orn = [];
for i = length(od):-1:1
	SNSdata.LFP = [SNSdata.LFP  od(i).LFP];
	SNSdata.orn = [SNSdata.orn; vectorise(i*ones(od(i).n_trials,1))];
	SNSdata.fA = [SNSdata.fA  od(i).firing_rate];
	SNSdata.PID = [SNSdata.PID  od(i).stimulus];
end

clear od
load('/local-data/DA-paper/data-for-paper/fig7/nat-stim-ab3/combined_data.ORNData','-mat')
DNSdata.LFP = [];
DNSdata.fA = [];
DNSdata.PID = [];
DNSdata.orn = [];
for i = length(od):-1:1
	DNSdata.LFP = [DNSdata.LFP  od(i).LFP];
	DNSdata.orn = [DNSdata.orn; vectorise(i*ones(od(i).n_trials,1))];
	DNSdata.fA = [DNSdata.fA  od(i).firing_rate];
	DNSdata.PID = [DNSdata.PID  od(i).stimulus];
end
% remove baseline from stimulus
for i = 1:size(DNSdata.PID,2)
	DNSdata.PID(:,i) = DNSdata.PID(:,i) - min(DNSdata.PID(1:5e3,i));
end

 ;;;;;;  ;;;;;;;;     ;;;    ;;;;;;;;   ;;;;;;  ;;;;;;;;    ;;    ;;    ;;;    ;;;;;;;; 
;;    ;; ;;     ;;   ;; ;;   ;;     ;; ;;    ;; ;;          ;;;   ;;   ;; ;;      ;;    
;;       ;;     ;;  ;;   ;;  ;;     ;; ;;       ;;          ;;;;  ;;  ;;   ;;     ;;    
 ;;;;;;  ;;;;;;;;  ;;     ;; ;;;;;;;;   ;;;;;;  ;;;;;;      ;; ;; ;; ;;     ;;    ;;    
      ;; ;;        ;;;;;;;;; ;;   ;;         ;; ;;          ;;  ;;;; ;;;;;;;;;    ;;    
;;    ;; ;;        ;;     ;; ;;    ;;  ;;    ;; ;;          ;;   ;;; ;;     ;;    ;;    
 ;;;;;;  ;;        ;;     ;; ;;     ;;  ;;;;;;  ;;;;;;;;    ;;    ;; ;;     ;;    ;;    

;;;;;;;; ;;;; ;;;;;;;;  ;;;; ;;    ;;  ;;;;;;      ;;;;;;;;     ;;;    ;;;;;;;; ;;;;;;;; 
;;        ;;  ;;     ;;  ;;  ;;;   ;; ;;    ;;     ;;     ;;   ;; ;;      ;;    ;;       
;;        ;;  ;;     ;;  ;;  ;;;;  ;; ;;           ;;     ;;  ;;   ;;     ;;    ;;       
;;;;;;    ;;  ;;;;;;;;   ;;  ;; ;; ;; ;;   ;;;;    ;;;;;;;;  ;;     ;;    ;;    ;;;;;;   
;;        ;;  ;;   ;;    ;;  ;;  ;;;; ;;    ;;     ;;   ;;   ;;;;;;;;;    ;;    ;;       
;;        ;;  ;;    ;;   ;;  ;;   ;;; ;;    ;;     ;;    ;;  ;;     ;;    ;;    ;;       
;;       ;;;; ;;     ;; ;;;; ;;    ;;  ;;;;;;      ;;     ;; ;;     ;;    ;;    ;;;;;;;; 


%% Firing rate: Sparse naturalistic stimulus
% First, I fit a NLN model to the firing rate responses to sparse naturalistic stimuli. The following figure shows the best fit model responses, and I also plot the 

clear p
p.   k0 = 0.1386;
p.tau_z = 10;
p.    B = 0;
p.  n_z = 1;
p.    n = 1;
p. tau1 = 24.2813;
p. tau2 = 78.3750;
p.  n_y = 2;
p.    A = 0.5356;
p.    C = 176.2152;

S = mean(SNSdata.PID(:,SNSdata.orn == 1),2);
R = mean(SNSdata.fA(:,SNSdata.orn == 1),2);
Rhat = aNLN2(S,p);
t = 1e-3*(1:length(R));

figure('outerposition',[0 0 1409 801],'PaperUnits','points','PaperSize',[1409 801]); hold on
subplot(2,4,1:4); hold on
plot(t,R,'k')
plot(t,Rhat,'r')
xlabel('Time (s)')
ylabel('ab3A Firing rate (Hz)')
legend({'ab3A','NLN model'})

subplot(2,4,5); hold on
x = logspace(-4,log10(max(S)*10),100);
a = hill([1 p.k0 p.n],x);
plot(x/mean(S),a,'r')
xlabel('S/<S>')
ylabel('a')
set(gca,'XScale','log','XLim',[1e-3 1e2],'XTick',[1e-3  1e-2 1e-1 1 10])

subplot(2,4,6); hold on
t = 1:600;
clear q
q.A = p.A; q.tau1 = p.tau1; q.tau2 = p.tau2; q.n = p.n_y;
K = filter_gamma2(t,q);
plot(t,K,'r')
xlabel('Filter lag (ms)')
ylabel('Filter')

subplot(2,4,7); hold on
plot(Rhat,R,'k.')
xlabel('NLN prediction (Hz)')
ylabel('ab3A firing rate (Hz)')
legend(['r^2 = ' oval(rsquare(R,Rhat))],'Location','northwest')

% vary k_D and n
all_k_D = logspace(log10(mean(S)/10),log10(mean(S)*10),31);
all_n = 1:10;
r2 = NaN(length(all_k_D),length(all_n));
if exist('.cache/sparse_nat_stim_k_D_n.mat','file')
	load('.cache/sparse_nat_stim_k_D_n.mat')
else
	for i = 1:length(all_k_D)
		textbar(i,length(all_k_D))
		for j = 1:length(all_n)
			q = p;
			q.k0 = all_k_D(i);
			q.n = all_n(j);
			Rhat = aNLN2(S,q);
			r2(i,j) = rsquare(Rhat,R);
		end
	end
	save('.cache/sparse_nat_stim_k_D_n.mat','r2')
end


n_labels = {};
for i = 1:2:length(all_n)
	n_labels{i} = oval(all_n(i));
end
k_D_labels = {'1/10','1','10'};

subplot(2,4,8); hold on
imagesc(r2)
ylabel('K_D/<S>')
xlabel('n')
colorbar
caxis([0.5 1])
set(gca,'XTick',[1:2:10],'XTickLabel',n_labels(cellfun(@(x) ~isempty(x),(n_labels))),'XTickLabelRotation',45)
set(gca,'YTick',[1 16 31],'YTickLabel',k_D_labels)
set(gca,'XLim',[.6 10.5],'YLim',[.5 31.5])
title('r^2(data, NLN model)')


prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%% 
% Does the addition of a adapting front-end improve the fit? I now allow the $k_D$ of the input nonlinearity to change with the stimulus over some recent history. 

clear p
p.   k0 = 0.1210;
p.tau_z = 8.0586;
p.    B = 0.0664;
p.  n_z = 9.8750;
p.    n = 1;
p. tau1 = 26.4727;
p. tau2 = 83.9062;
p.  n_y = 2;
p.    A = 0.5513;
p.    C = 187.5902;

figure('outerposition',[0 0 1000 803],'PaperUnits','points','PaperSize',[1000 803]); hold on

% solve
[Rhat,Ky,Kz,k_D] = aNLN2(S,p);

subplot(2,3,1); hold on
all_k_D = [min(k_D) mean(k_D) max(k_D)];
c = parula(3);
for i = 1:length(all_k_D)
	x = logspace(-4,log10(max(S)*10),100);
	a = hill([1 all_k_D(i) p.n],x);
	plot(x/mean(S),a,'Color',c(i,:))
end
xlabel('S/<S>')
ylabel('a')
set(gca,'XScale','log','XLim',[1e-3 1e2],'XTick',[1e-3  1e-2 1e-1 1 10])
legend({'min k_D','mean k_D','max k_D'},'Location','northwest')

subplot(2,3,2); hold on
plot(1:length(Ky),Ky,'r')
plot(1:length(Kz),Kz,'b')
legend({'Response filter','Adaptation filter'})
xlabel('Lag (ms)')

subplot(2,3,3); hold on
plot(Rhat,R,'k.')
xlabel('NLN prediction (Hz)')
ylabel('ab3A firing rate (Hz)')
legend(['r^2 = ' oval(rsquare(R,Rhat))],'Location','northwest')

subplot(2,3,4); hold on
[hy,hx] = histcounts(k_D,100);
hy = hy/sum(hy);
plot(hx(2:end)/mean(S),hy,'k')
set(gca,'XScale','log','YScale','log')
xlabel('k_D/<S>')
ylabel('Probability')

% vary B and tau_gain
all_tau_z = unique(round(logspace(0,3,31)));
all_B = unique(logspace(-3,3,31));
r2 = NaN(length(all_tau_z),length(all_n));
if exist('.cache/sparse_nat_stim_B_tau.mat','file')
	load('.cache/sparse_nat_stim_B_tau.mat')
else
	for i = 1:length(all_tau_z)
		textbar(i,length(all_tau_z))
		for j = 1:length(all_B)
			q = p;
			q.tau_z = all_tau_z(i);
			q.B = all_B(j);
			Rhat = aNLN2(S,q);
			r2(i,j) = rsquare(Rhat,R);
		end
	end
	save('.cache/sparse_nat_stim_B_tau.mat','r2')
end


B_labels = {};
b_tick = false(length(all_B),1);
for i = 1:length(all_B)
	if log10(all_B(i)) == round(log10(all_B(i)))
		B_labels{i} = ['10^{' oval(log10(all_B(i))) '}'];
		b_tick(i) = true;
	end
end


tau_labels = {};
tau_tick = false(length(all_tau_z),1);
for i = 1:length(all_tau_z)
	if log10(all_tau_z(i)) == round(log10(all_tau_z(i)))
		tau_labels{i} = ['10^{' oval(log10(all_tau_z(i))) '}'];
		tau_tick(i) = true;
	end
end

subplot(2,3,5:6); hold on
imagesc(r2)
ylabel('\tau_{adaptation} (ms)')
xlabel('\beta')
colorbar
caxis([0.5 1])
set(gca,'XTick',find(b_tick),'XTickLabel',B_labels(b_tick),'XTickLabelRotation',45)
set(gca,'YTick',find(tau_tick),'YTickLabel',tau_labels(tau_tick),'YTickLabelRotation',45)
set(gca,'XLim',[.5 size(r2,2)+.5],'YLim',[.5 size(r2,1)+.5])
title('r^2(data, NLN model)')

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

;;;;;;;;  ;;;;;;;; ;;    ;;  ;;;;;;  ;;;;;;;;    ;;    ;;    ;;;    ;;;;;;;; 
;;     ;; ;;       ;;;   ;; ;;    ;; ;;          ;;;   ;;   ;; ;;      ;;    
;;     ;; ;;       ;;;;  ;; ;;       ;;          ;;;;  ;;  ;;   ;;     ;;    
;;     ;; ;;;;;;   ;; ;; ;;  ;;;;;;  ;;;;;;      ;; ;; ;; ;;     ;;    ;;    
;;     ;; ;;       ;;  ;;;;       ;; ;;          ;;  ;;;; ;;;;;;;;;    ;;    
;;     ;; ;;       ;;   ;;; ;;    ;; ;;          ;;   ;;; ;;     ;;    ;;    
;;;;;;;;  ;;;;;;;; ;;    ;;  ;;;;;;  ;;;;;;;;    ;;    ;; ;;     ;;    ;;    

;;;;;;;; ;;;; ;;;;;;;;  ;;;; ;;    ;;  ;;;;;;      ;;;;;;;;     ;;;    ;;;;;;;; ;;;;;;;; 
;;        ;;  ;;     ;;  ;;  ;;;   ;; ;;    ;;     ;;     ;;   ;; ;;      ;;    ;;       
;;        ;;  ;;     ;;  ;;  ;;;;  ;; ;;           ;;     ;;  ;;   ;;     ;;    ;;       
;;;;;;    ;;  ;;;;;;;;   ;;  ;; ;; ;; ;;   ;;;;    ;;;;;;;;  ;;     ;;    ;;    ;;;;;;   
;;        ;;  ;;   ;;    ;;  ;;  ;;;; ;;    ;;     ;;   ;;   ;;;;;;;;;    ;;    ;;       
;;        ;;  ;;    ;;   ;;  ;;   ;;; ;;    ;;     ;;    ;;  ;;     ;;    ;;    ;;       
;;       ;;;; ;;     ;; ;;;; ;;    ;;  ;;;;;;      ;;     ;; ;;     ;;    ;;    ;;;;;;;; 


%% Firing rate: dense naturalistic stimulus
% Now I fit a non-adapting NLN model to the dense naturalistic stimulus. 

S = nanmean(DNSdata.PID(:,DNSdata.orn == 6),2);
R = nanmean(DNSdata.fA(:,DNSdata.orn == 6),2);

clear p
p.   k0 = 0.0656;
p.tau_z = 1;
p.    B = 0;
p.  n_z = 1;
p.    n = 1.3833;
p. tau1 = 23.3110;
p. tau2 = 44.3750;
p.  n_y = 2;
p.    A = 0.7908;
p.    C = 316.1258;


Rhat = aNLN2(S,p);
t = 1e-3*(1:length(R));

figure('outerposition',[0 0 1409 801],'PaperUnits','points','PaperSize',[1409 801]); hold on
subplot(2,4,1:4); hold on
plot(t,R,'k')
plot(t,Rhat,'r')
xlabel('Time (s)')
ylabel('ab3A Firing rate (Hz)')
legend({'ab3A','NLN model'})

subplot(2,4,5); hold on
x = logspace(-4,log10(max(S)*10),100);
a = hill([1 p.k0 p.n],x);
plot(x/mean(S),a,'r')
xlabel('S/<S>')
ylabel('a')
set(gca,'XScale','log','XLim',[1e-3 1e2],'XTick',[1e-3  1e-2 1e-1 1 10])

subplot(2,4,6); hold on
t = 1:600;
clear q
q.A = p.A; q.tau1 = p.tau1; q.tau2 = p.tau2; q.n = p.n_y;
K = filter_gamma2(t,q);
plot(t,K,'r')
xlabel('Filter lag (ms)')
ylabel('Filter')

subplot(2,4,7); hold on
plot(Rhat,R,'k.')
xlabel('NLN prediction (Hz)')
ylabel('ab3A firing rate (Hz)')
legend(['r^2 = ' oval(rsquare(R,Rhat))],'Location','northwest')

% vary k_D and n
all_k_D = logspace(log10(mean(S)/10),log10(mean(S)*10),31);
all_n = 1:10;
r2 = NaN(length(all_k_D),length(all_n));
if exist('.cache/dense_nat_stim_k_D_n.mat','file')
	load('.cache/dense_nat_stim_k_D_n.mat')
else
	for i = 1:length(all_k_D)
		textbar(i,length(all_k_D))
		for j = 1:length(all_n)
			q = p;
			q.k0 = all_k_D(i);
			q.n = all_n(j);
			Rhat = aNLN2(S,q);
			r2(i,j) = rsquare(Rhat,R);
		end
	end
	save('.cache/dense_nat_stim_k_D_n.mat','r2')
end


n_labels = {};
for i = 1:2:length(all_n)
	n_labels{i} = oval(all_n(i));
end
k_D_labels = {'1/10','1','10'};

subplot(2,4,8); hold on
imagesc(r2)
ylabel('K_D/<S>')
xlabel('n')
colorbar
caxis([0 1])
set(gca,'XTick',[1:2:10],'XTickLabel',n_labels(cellfun(@(x) ~isempty(x),(n_labels))),'XTickLabelRotation',45)
set(gca,'YTick',[1 16 31],'YTickLabel',k_D_labels)
set(gca,'XLim',[.6 10.5],'YLim',[.5 31.5])
title('r^2(data, NLN model)')


prettyFig();

if being_published
	snapnow
	delete(gcf)
end


%% 
% Does the addition of a adapting front-end improve the fit? I now allow the $k_D$ of the input nonlinearity to change with the stimulus over some recent history. 

clear p
p.   k0 = 0.0656;
p.tau_z = 1000;
p.    B = 0.0391;
p.  n_z = 2;
p.    n = 1.3799;
p. tau1 = 23.6318;
p. tau2 = 41.0430;
p.  n_y = 2;
p.    A = 0.7906;
p.    C = 345.0213;

figure('outerposition',[0 0 1000 803],'PaperUnits','points','PaperSize',[1000 803]); hold on

% solve
[Rhat,Ky,Kz,k_D] = aNLN2(S,p);

subplot(2,3,1); hold on
all_k_D = [min(k_D) mean(k_D) max(k_D)];
c = parula(3);
for i = 1:length(all_k_D)
	x = logspace(-4,log10(max(S)*10),100);
	a = hill([1 all_k_D(i) p.n],x);
	plot(x/mean(S),a,'Color',c(i,:))
end
xlabel('S/<S>')
ylabel('a')
set(gca,'XScale','log','XLim',[1e-3 1e2],'XTick',[1e-3  1e-2 1e-1 1 10])
legend({'min k_D','mean k_D','max k_D'},'Location','northwest')

subplot(2,3,2); hold on
plot(1:length(Ky),Ky,'r')
plot(1:length(Kz),Kz,'b')
legend({'Response filter','Adaptation filter'})
xlabel('Lag (ms)')

subplot(2,3,3); hold on
plot(Rhat,R,'k.')
xlabel('NLN prediction (Hz)')
ylabel('ab3A firing rate (Hz)')
legend(['r^2 = ' oval(rsquare(R,Rhat))],'Location','northwest')

subplot(2,3,4); hold on
[hy,hx] = histcounts(k_D,100);
hy = hy/sum(hy);
plot(hx(2:end)/mean(S),hy,'k')
set(gca,'XScale','log','YScale','linear','XLim',[.1 10])
xlabel('k_D/<S>')
ylabel('Probability')

% vary B and tau_gain
all_tau_z = unique(round(logspace(0,3,31)));
all_B = unique(logspace(-3,3,31));
r2 = NaN(length(all_tau_z),length(all_n));
if exist('.cache/dense_nat_stim_B_tau.mat','file')
	load('.cache/dense_nat_stim_B_tau.mat')
else
	for i = 1:length(all_tau_z)
		textbar(i,length(all_tau_z))
		for j = 1:length(all_B)
			q = p;
			q.tau_z = all_tau_z(i);
			q.B = all_B(j);
			Rhat = aNLN2(S,q);
			r2(i,j) = rsquare(Rhat,R);
		end
	end
	save('.cache/dense_nat_stim_B_tau.mat','r2')
end


B_labels = {};
b_tick = false(length(all_B),1);
for i = 1:length(all_B)
	if log10(all_B(i)) == round(log10(all_B(i)))
		B_labels{i} = ['10^{' oval(log10(all_B(i))) '}'];
		b_tick(i) = true;
	end
end


tau_labels = {};
tau_tick = false(length(all_tau_z),1);
for i = 1:length(all_tau_z)
	if log10(all_tau_z(i)) == round(log10(all_tau_z(i)))
		tau_labels{i} = ['10^{' oval(log10(all_tau_z(i))) '}'];
		tau_tick(i) = true;
	end
end

subplot(2,3,5:6); hold on
imagesc(r2)
ylabel('\tau_{adaptation} (ms)')
xlabel('\beta')
colorbar
caxis([0 1])
set(gca,'XTick',find(b_tick),'XTickLabel',B_labels(b_tick),'XTickLabelRotation',45)
set(gca,'YTick',find(tau_tick),'YTickLabel',tau_labels(tau_tick),'YTickLabelRotation',45)
set(gca,'XLim',[.5 size(r2,2)+.5],'YLim',[.5 size(r2,1)+.5])
title('r^2(data, NLN model)')

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

;;     ;;  ;;;;;;   ;;;;;;   
;;;   ;;; ;;    ;; ;;    ;;  
;;;; ;;;; ;;       ;;        
;; ;;; ;;  ;;;;;;  ;;   ;;;; 
;;     ;;       ;; ;;    ;;  
;;     ;; ;;    ;; ;;    ;;  
;;     ;;  ;;;;;;   ;;;;;;   

;;;;;;;; ;;;; ;;;;;;;;  ;;;; ;;    ;;  ;;;;;;      ;;;;;;;;     ;;;    ;;;;;;;; ;;;;;;;; 
;;        ;;  ;;     ;;  ;;  ;;;   ;; ;;    ;;     ;;     ;;   ;; ;;      ;;    ;;       
;;        ;;  ;;     ;;  ;;  ;;;;  ;; ;;           ;;     ;;  ;;   ;;     ;;    ;;       
;;;;;;    ;;  ;;;;;;;;   ;;  ;; ;; ;; ;;   ;;;;    ;;;;;;;;  ;;     ;;    ;;    ;;;;;;   
;;        ;;  ;;   ;;    ;;  ;;  ;;;; ;;    ;;     ;;   ;;   ;;;;;;;;;    ;;    ;;       
;;        ;;  ;;    ;;   ;;  ;;   ;;; ;;    ;;     ;;    ;;  ;;     ;;    ;;    ;;       
;;       ;;;; ;;     ;; ;;;; ;;    ;;  ;;;;;;      ;;     ;; ;;     ;;    ;;    ;;;;;;;; 


%% Firing rate: mean shifted Gaussians
% Now, I attempt to fit a non-adapting NLN model to the mean shifted gaussian data. 

c = 1;
for i = [2:10]
	R = MSGdata.fA(35e3:55e3,MSGdata.paradigm == i);
	S = MSGdata.PID(35e3:55e3,MSGdata.paradigm == i);
	S(:,sum(R)==0) = [];
	R(:,sum(R)==0) = [];
	data(c).response = nanmean(R,2);
	data(c).response(1:2e3) = NaN;
	data(c).stimulus = nanmean(S,2);
	c = c + 1;
end

clear p
p.   k0 = 0.2994;
p.tau_z = 1;
p.    B = 0;
p.  n_z = 1;
p.    n = 3.8003;
p. tau1 = 32.6495;
p. tau2 = 40.0469;
p.  n_y = 2;
p.    A = 0.7906;
p.    C = 225.1171;

% generate responses 
for i = 1:size(MSGdata.PID,2)
	MSGdata.NLN_R(:,i) = aNLN2(MSGdata.PID(:,i),p);
end

figure('outerposition',[0 0 1300 901],'PaperUnits','points','PaperSize',[1300 901]); hold on
c = parula(11);

time = 1e-3*(1:length(MSGdata.PID));

show_these_paradigms = [2 6 8 10];
for i = 1:4
	subplot(2,4,i); hold on
	R = MSGdata.fA(:,MSGdata.paradigm == show_these_paradigms(i));
	X = MSGdata.NLN_R(:,MSGdata.paradigm == show_these_paradigms(i));
	R(:,sum(R)==0) = [];
	R = nanmean(R,2);
	X = nanmean(X,2);
	plot(time,R,'k')
	plot(time,X,'r')
	set(gca,'XLim',[40 50])
end

subplot(2,4,5); hold on
x = logspace(-2,1,100);
a = hill([1 p.k0 p.n],x);
plot(x,a,'r')
xlabel('S')
ylabel('a')
set(gca,'XScale','log','XLim',[1e-2 1e1],'XTick',[1e-3  1e-2 1e-1 1 10])

[~,Ky] = aNLN2(MSGdata.PID(:,1),p);
subplot(2,4,6); hold on
filtertime = 1e-3*(1:length(Ky));
plot(filtertime,Ky,'r')
xlabel('Lag (s)')
ylabel('Filter')

subplot(2,4,7); hold on
all_X = []; all_R = [];
for i = 1:10
	R = MSGdata.fA(35e3:55e3,MSGdata.paradigm == i);
	X = MSGdata.NLN_R(35e3:55e3,MSGdata.paradigm == i);
	R(:,sum(R)==0) = [];
	X(:,sum(X)==0) = [];
	R = nanmean(R,2);
	X = nanmean(X,2);
	plot(X(1:50:end),R(1:50:end),'.','Color',c(i,:));
	all_X = [all_X; X]; all_R = [all_R; R];
end


% vary k_D and n
all_k_D = logspace(log10(p.k0/10),log10(p.k0*10),31);
all_n = 1:10;
r2 = NaN(length(all_k_D),length(all_n));
if exist('.cache/msg_k_D_n.mat','file')
	load('.cache/msg_k_D_n.mat')
else
	for i = 1:length(all_k_D)
		textbar(i,length(all_k_D))
		for j = 1:length(all_n)
			all_X = []; 
			Rhat = MSGdata.fA;
			q = p;
			q.k0 = all_k_D(i);
			q.n = all_n(j);
			for k = 1:size(MSGdata.PID,2)
				Rhat(:,k) = aNLN2(MSGdata.PID(:,k),q);
			end

			for k = 1:10
				R = MSGdata.fA(35e3:55e3,MSGdata.paradigm == k);
				X = Rhat(35e3:55e3,MSGdata.paradigm == k);
				R(:,sum(R)==0) = [];
				X(:,sum(X)==0) = [];
				R = nanmean(R,2);
				X = nanmean(X,2);
				all_X = [all_X; X]; 
			end
			r2(i,j) = rsquare(all_R, all_X);
		end
	end
	save('.cache/msg_k_D_n.mat','r2')
end


n_labels = {};
for i = 1:2:length(all_n)
	n_labels{i} = oval(all_n(i));
end
k_D_labels = {'1/10','1','10'};

subplot(2,4,8); hold on
imagesc(r2)
ylabel('K_D/k_{D best fit}')
xlabel('n')
colorbar
caxis([0 1])
set(gca,'XTick',[1:2:10],'XTickLabel',n_labels(cellfun(@(x) ~isempty(x),(n_labels))),'XTickLabelRotation',45)
set(gca,'YTick',[1 16 31],'YTickLabel',k_D_labels)
set(gca,'XLim',[.6 10.5],'YLim',[.5 31.5])
title('r^2(data, NLN model)')


prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% Now I fit an adapting NLN model to the same data, and see if we can capture the observed change in gain. 

clear p
p.   k0 = 0.1080;
p.tau_z = 44.3828;
p.    B = 0.3164;
p.  n_z = 1.0312;
p.    n = 2.2026;
p. tau1 = 32.0323;
p. tau2 = 199.9961;
p.  n_y = 2;
p.    A = 0.7906;
p.    C = 144.6366;


% generate responses 
for i = 1:size(MSGdata.PID,2)
	MSGdata.NLN_R(:,i) = aNLN2(MSGdata.PID(:,i),p);
end

figure('outerposition',[0 0 1400 901],'PaperUnits','points','PaperSize',[1400 901]); hold on
c = parula(11);

time = 1e-3*(1:length(MSGdata.PID));

show_these_paradigms = [2 4 6 8];
for i = 1:4
	subplot(2,4,i); hold on
	R = MSGdata.fA(:,MSGdata.paradigm == show_these_paradigms(i));
	X = MSGdata.NLN_R(:,MSGdata.paradigm == show_these_paradigms(i));
	R(:,sum(R)==0) = [];
	R = nanmean(R,2);
	X = nanmean(X,2);
	plot(time,R,'k')
	plot(time,X,'r')
	set(gca,'XLim',[45 55],'YLim',[0 70])
	xlabel('Time (s)')
	ylabel('Firing rate (Hz)')
end

subplot(2,4,5); hold on
x = logspace(-2,1,100);
for i = 1:10
	mean_k_D = [];
	for j = 1:size(MSGdata.PID,2)
		if MSGdata.paradigm(j) == i
			[~,~,~,k_D] = aNLN2(MSGdata.PID(:,j),p);
			mean_k_D = [mean_k_D k_D(35e3:55e3)];
		end
	end
	mean_k_D = mean(mean(mean_k_D));
	a = hill([1 mean_k_D p.n],x);
	plot(x,a,'Color',c(i,:))
end
xlabel('S')
ylabel('a')
set(gca,'XScale','log','XLim',[1e-2 1e1],'XTick',[1e-3  1e-2 1e-1 1 10])

[~,Ky,Kz] = aNLN2(MSGdata.PID(:,1),p);
subplot(2,4,6); hold on
plot(1:length(Ky),Ky,'r')
plot(1:length(Kz),Kz,'b')
legend({'Response filter','Adaptation filter'})
xlabel('Lag (ms)')
ylabel('Filter')


subplot(2,4,7); hold on
all_X = []; all_R = [];
for i = 1:10
	R = MSGdata.fA(35e3:55e3,MSGdata.paradigm == i);
	X = MSGdata.NLN_R(35e3:55e3,MSGdata.paradigm == i);
	R(:,sum(R)==0) = [];
	X(:,sum(X)==0) = [];
	R = nanmean(R,2);
	X = nanmean(X,2);
	plot(X(1:50:end),R(1:50:end),'.','Color',c(i,:));
	all_X = [all_X; X]; all_R = [all_R; R];
end
xlabel('model prediction (Hz)')
ylabel('ab3A firing rate (Hz)')


% vary tau_gain and B
all_tau_z = unique(round(logspace(0,4,21)));
all_B = unique(logspace(-3,3,31));
r2 = NaN(length(all_tau_z),length(all_B));
if exist('.cache/msg_tau_B.mat','file')
	load('.cache/msg_tau_B.mat')
else
	for i = 1:length(all_tau_z)
		textbar(i,length(all_tau_z))
		for j = 1:length(all_B)
			all_X = []; 
			Rhat = MSGdata.fA;
			q = p;
			q.tau_z = all_tau_z(i);
			q.B = all_B(j);
			for k = 1:size(MSGdata.PID,2)
				Rhat(:,k) = aNLN2(MSGdata.PID(:,k),q);
			end

			for k = 1:10
				R = MSGdata.fA(35e3:55e3,MSGdata.paradigm == k);
				X = Rhat(35e3:55e3,MSGdata.paradigm == k);
				R(:,sum(R)==0) = [];
				X(:,sum(R)==0) = [];
				R = nanmean(R,2);
				X = nanmean(X,2);
				all_X = [all_X; X]; 
			end
			r2(i,j) = rsquare(all_R, all_X);
		end
	end
	save('.cache/msg_tau_B.mat','r2')
end

B_labels = {};
b_tick = false(length(all_B),1);
for i = 1:length(all_B)
	if log10(all_B(i)) == round(log10(all_B(i)))
		B_labels{i} = ['10^{' oval(log10(all_B(i))) '}'];
		b_tick(i) = true;
	end
end


tau_labels = {};
tau_tick = false(length(all_tau_z),1);
for i = 1:length(all_tau_z)
	if log10(all_tau_z(i)) == round(log10(all_tau_z(i)))
		tau_labels{i} = ['10^{' oval(log10(all_tau_z(i))) '}'];
		tau_tick(i) = true;
	end
end

ax = subplot(2,4,8); hold on
imagesc(r2)
ylabel('\tau_{adaptation} (ms)')
xlabel('\beta')
colorbar
caxis([0 1])
set(gca,'XTick',find(b_tick),'XTickLabel',B_labels(b_tick),'XTickLabelRotation',45)
set(gca,'YTick',find(tau_tick),'YTickLabel',tau_labels(tau_tick),'YTickLabelRotation',45)
set(gca,'XLim',[.5 size(r2,2)+.5],'YLim',[.5 size(r2,1)+.5])
title('r^2(data, NLN model)')
ax.Position = [0.7619    0.1100    0.15    0.3412];

prettyFig();

if being_published
	snapnow
	delete(gcf)
end


%% Firing rate: variance gain control
% 

%% Firing rate: All data
% 




%% Version Info
%
pFooter;


