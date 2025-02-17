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
load(getPath(dataManager,'0f4e4e0686913ee3c4fa450d9e2bb344'),'-mat')
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
load(getPath(dataManager,'aeb361c027b71938021c12a6a12a85cd'),'-mat')
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

% get the variance switching data
clear VSdata
[VSdata.PID, VSdata.LFP, VSdata.fA, VSdata.paradigm, VSdata.orn, VSdata.fly] = consolidateData(getPath(dataManager,'e30707e8e8ef6c0d832eee31eaa585aa'),1);
% remove baseline from stimulus
VSdata.PID = bsxfun(@minus, VSdata.PID, min(VSdata.PID));
VSdata.LFP = bsxfun(@minus, VSdata.LFP, mean(VSdata.LFP(1:4000,:)));

% get the variance switching data filters 
load(getPath(dataManager,'457ee16a326f47992e35a7d5281f9cc4'));
VSdata.K1 = nanmean(K1,2);

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
% First, I fit a NLN model to the firing rate responses to sparse naturalistic stimuli. The following figure shows the ORN firing rate with the best-fit NLN model (top). In the next two panels, I plot the best-fit input nonlinearity and the best-fit filter (together, the NL model). I compare the data to the best-fit prediction, and show the $r^2$ of the prediction (it's very good). Finally, I vary two parameters of the input nonlinearity (it is a Hill function), to determine how strongly this data constrains the shape of the input nonlinearity (it looks like it doesn't constrain it very tightly). 

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

                         ;;;    ;;;;;;;;     ;;;    ;;;;;;;;  ;;;;;;;; 
           ;;           ;; ;;   ;;     ;;   ;; ;;   ;;     ;;    ;;    
           ;;          ;;   ;;  ;;     ;;  ;;   ;;  ;;     ;;    ;;    
         ;;;;;;       ;;     ;; ;;     ;; ;;     ;; ;;;;;;;;     ;;    
           ;;         ;;;;;;;;; ;;     ;; ;;;;;;;;; ;;           ;;    
           ;;         ;;     ;; ;;     ;; ;;     ;; ;;           ;;    
                      ;;     ;; ;;;;;;;;  ;;     ;; ;;           ;;    




%% 
% Does the addition of a adapting front-end improve the fit? I now allow the $k_D$ of the input nonlinearity to change with the stimulus over some recent history, and then fit this model to the same data. First, I plot the extememum and mean values of $k_D$ from this model, to show how the effect of this adaptation on the model (top left). I then show the response and "adaptation" filters, and then compare the data to the model prediction (top right). I also plot the probability distribution of the dynamically-updating $k_D$, and we see that it barely changes from its mean value. Finally, I vary the parameters of update of $k_D$, the timescale and the degree, and show that this data does not strongly constrain the timescale of front-end gain control (nor does it really need it). 

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
xlabel('adapting NLN prediction (Hz)')
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
caxis([0 1])
set(gca,'XTick',find(b_tick),'XTickLabel',B_labels(b_tick),'XTickLabelRotation',45)
set(gca,'YTick',find(tau_tick),'YTickLabel',tau_labels(tau_tick),'YTickLabelRotation',45)
set(gca,'XLim',[.5 size(r2,2)+.5],'YLim',[.5 size(r2,1)+.5])
title('r^2(data, adapting NLN model)')

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

                         ;;;    ;;;;;;;;     ;;;    ;;;;;;;;  ;;;;;;;; 
           ;;           ;; ;;   ;;     ;;   ;; ;;   ;;     ;;    ;;    
           ;;          ;;   ;;  ;;     ;;  ;;   ;;  ;;     ;;    ;;    
         ;;;;;;       ;;     ;; ;;     ;; ;;     ;; ;;;;;;;;     ;;    
           ;;         ;;;;;;;;; ;;     ;; ;;;;;;;;; ;;           ;;    
           ;;         ;;     ;; ;;     ;; ;;     ;; ;;           ;;    
                      ;;     ;; ;;;;;;;;  ;;     ;; ;;           ;;    




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
xlabel('adapting NLN prediction (Hz)')
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
title('r^2(data, adapting NLN model)')

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
% Now, I attempt to fit a non-adapting NLN model to the mean shifted gaussian data. In the top row, I plot the responses of the neuron in black to four different mean stimuli, and compare it to the model predictions in red. We clearly see that 1) the model tends to over-estimate mean firing rates as the mean stimulus increases, a consequence of saturation, and 2) the fluctuations in the model response fall of much faster with increasing mean stimulus than in the data. This is also visible in the comparison of ORN responses to model responses for the different mean stimuli. 

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

figure('outerposition',[0 0 1300 744],'PaperUnits','points','PaperSize',[1300 744]); hold on
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
	set(gca,'XLim',[40 50],'YLim',[0 70])
	if i == 1
		xlabel('Time (s)')
		ylabel('Firing rate (Hz)')
	end
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
xlabel('Model prediction (Hz)')
ylabel('ab3A firing rate (Hz)')

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
      
                         ;;;    ;;;;;;;;     ;;;    ;;;;;;;;  ;;;;;;;; 
           ;;           ;; ;;   ;;     ;;   ;; ;;   ;;     ;;    ;;    
           ;;          ;;   ;;  ;;     ;;  ;;   ;;  ;;     ;;    ;;    
         ;;;;;;       ;;     ;; ;;     ;; ;;     ;; ;;;;;;;;     ;;    
           ;;         ;;;;;;;;; ;;     ;; ;;;;;;;;; ;;           ;;    
           ;;         ;;     ;; ;;     ;; ;;     ;; ;;           ;;    
                      ;;     ;; ;;;;;;;;  ;;     ;; ;;           ;;    



%%
% Now I fit an adapting NLN model to the same data, and see if we can capture the observed change in gain. We see that this modification allows the model to change its input nonlinearity, moving it to the right with increasing mean stimulus, and allows it to approximate the ORN response much better. However, we see that this data poorly constrains the timescale of gain control in this model. 

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
title('r^2(data, adapting NLN model)')
ax.Position = [0.7619    0.1100    0.15    0.3412];

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

;;     ;;    ;;;    ;;;;;;;;  ;;;;    ;;;    ;;    ;;  ;;;;;;  ;;;;;;;; 
;;     ;;   ;; ;;   ;;     ;;  ;;    ;; ;;   ;;;   ;; ;;    ;; ;;       
;;     ;;  ;;   ;;  ;;     ;;  ;;   ;;   ;;  ;;;;  ;; ;;       ;;       
;;     ;; ;;     ;; ;;;;;;;;   ;;  ;;     ;; ;; ;; ;; ;;       ;;;;;;   
 ;;   ;;  ;;;;;;;;; ;;   ;;    ;;  ;;;;;;;;; ;;  ;;;; ;;       ;;       
  ;; ;;   ;;     ;; ;;    ;;   ;;  ;;     ;; ;;   ;;; ;;    ;; ;;       
   ;;;    ;;     ;; ;;     ;; ;;;; ;;     ;; ;;    ;;  ;;;;;;  ;;;;;;;; 

;;;;;;;; ;;;; ;;;;;;;;  ;;;; ;;    ;;  ;;;;;;      ;;;;;;;;     ;;;    ;;;;;;;; ;;;;;;;; 
;;        ;;  ;;     ;;  ;;  ;;;   ;; ;;    ;;     ;;     ;;   ;; ;;      ;;    ;;       
;;        ;;  ;;     ;;  ;;  ;;;;  ;; ;;           ;;     ;;  ;;   ;;     ;;    ;;       
;;;;;;    ;;  ;;;;;;;;   ;;  ;; ;; ;; ;;   ;;;;    ;;;;;;;;  ;;     ;;    ;;    ;;;;;;   
;;        ;;  ;;   ;;    ;;  ;;  ;;;; ;;    ;;     ;;   ;;   ;;;;;;;;;    ;;    ;;       
;;        ;;  ;;    ;;   ;;  ;;   ;;; ;;    ;;     ;;    ;;  ;;     ;;    ;;    ;;       
;;       ;;;; ;;     ;; ;;;; ;;    ;;  ;;;;;;      ;;     ;; ;;     ;;    ;;    ;;;;;;;; 




%% Firing rate: variance gain control
% In this section, I fit a non-adapting NLN model to the variance switch data. The top panel shows the neuron's response with the best fit model. In the bottom left panel, I plot the neuron's firing rate vs. the best-fit NLN model prediction, seperately for the low-variance epochs (blue) and the high-variance epochs (red). Note that both curves have the same slope, indicating that the model has captured whatever gain change the neuron has done. This suggests that this data can precisely constrain some parameters of the input. In the bottom middle panel, I plot the ratio of gains (calcualted w.r.t to the NLN model) during the low and high variance epochs. When the curve crosses 1, the model shows the same degree of variance gain control that is seen in the data. I also plot the $r^2$ as a function of the steepness parameter. Note that the maximum of this plot does not occur at the same value of $n$ as when the model shows the same amount of variance gain control. 

example_trial = 4;
S = VSdata.PID(:,example_trial);
R = VSdata.fA(:,example_trial);

clear p
p.   k0 = 0.4115;
p.tau_z = 1;
p.    B = 0;
p.  n_z = 1;
p.    n = 3.9336;
p.n = 7;
p. tau1 = 14.1570;
p. tau2 = 44.6878;
p.  n_y = 4.6729;
p.    A = 0.3812;
p.    C = 86.2114;

X = aNLN2(S,p);
time = 1e-3*(1:length(X));

figure('outerposition',[0 0 1400 901],'PaperUnits','points','PaperSize',[1400 901]); hold on
subplot(2,1,1); hold on
plot(time,R,'k')
plot(time,X,'r')
set(gca,'XLim',[40 50])
xlabel('Time (s)')
ylabel('Firing rate (Hz)')

% reshape the response and the prediction
R = reshape(R,1e4,length(R)/1e4);
X = reshape(X,1e4,length(X)/1e4);
X(:,1) = []; X(:,end) = [];
R(:,1) = []; R(:,end) = [];

subplot(2,3,4); hold on
plotPieceWiseLinear(vectorise(X(1e3:5e3,:)),vectorise(R(1e3:5e3,:)),'nbins',40,'Color',[1 0 0]); 
plotPieceWiseLinear(vectorise(X(6e3:end,:)),vectorise(R(6e3:end,:)),'nbins',40,'Color',[0 0 1]); 
xlabel('NLN model prediction (Hz)')
ylabel('ab3A firing rate (Hz)')

% vary n, measure slopes and r^2
all_n = 1:21;
r2 = NaN(length(all_n),1);
slopes_lo = NaN(length(all_n),1);
slopes_hi = NaN(length(all_n),1);

if exist('.cache/variance_n.mat','file')
	load('.cache/variance_n.mat')
else
	for i = 1:length(all_n)
		textbar(i,length(all_n))
		q = p;
		q.n = all_n(i);

		Rhat = aNLN2(S,q);

		r2(i) = rsquare(VSdata.fA(1e4:end-1e4,example_trial), Rhat(1e4:end-1e4));

		Rhat = reshape(Rhat,1e4,length(Rhat)/1e4);
		Rhat(:,1) = []; Rhat(:,end) = [];

		ff = fit(vectorise(Rhat(1e3:5e3,:)),vectorise(R(1e3:5e3,:)),'poly1');
		slopes_hi(i) = ff.p1;

		ff = fit(vectorise(Rhat(6e3:end,:)),vectorise(R(6e3:end,:)),'poly1');
		slopes_lo(i) = ff.p1;
		
	end
	save('.cache/variance_n.mat','r2','slopes_lo','slopes_hi')
end

subplot(2,3,5); hold on
plot(all_n,slopes_lo(:)./slopes_hi(:),'k')
xlabel('n')
ylabel('gain_{low}/gain_{high}')
plot(all_n,0*all_n+1,'k--')

subplot(2,3,6); hold on
plot(all_n,r2,'k+-')
xlabel('n')
ylabel('r^2')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% Notice in the previous plot that even though the models got the gain in the high and low variance epochs right, it mis-estimated the mean response. This is because the mean stimulus changes slightly, and since there is no adaptation in this model, it doesn't get this right. What if I add in mean adaptation? We see that now it gets the slopes right, and the mean response right (top right). 

example_trial = 4;
S = VSdata.PID(:,example_trial);
R = VSdata.fA(:,example_trial);

clear p
p.    k0 = 1.0000e-06;
p. tau_z = 500;
p.B = .95;
p.   n_z = 1;
p.n = 7;
p.  tau1 = 15.8211;
p.  tau2 = 32.7425;
p.   n_y = 4.4971;
p.     A = 0.3851;
p.     C = 87.2973;
X = aNLN2(S,p);
time = 1e-3*(1:length(X));

figure('outerposition',[0 0 1450 802],'PaperUnits','points','PaperSize',[1450 802]); hold on
subplot(2,4,1:3); hold on
plot(time,R,'k')
plot(time,X,'r')
set(gca,'XLim',[40 50])
xlabel('Time (s)')
ylabel('Firing rate (Hz)')

% reshape the response and the prediction
R = reshape(R,1e4,length(R)/1e4);
X = reshape(X,1e4,length(X)/1e4);
X(:,1) = []; X(:,end) = [];
R(:,1) = []; R(:,end) = [];

subplot(2,4,4); hold on
plotPieceWiseLinear(vectorise(X(1e3:5e3,:)),vectorise(R(1e3:5e3,:)),'nbins',40,'Color',[1 0 0]); 
plotPieceWiseLinear(vectorise(X(6e3:end,:)),vectorise(R(6e3:end,:)),'nbins',40,'Color',[0 0 1]); 
xlabel('adapting NLN model prediction (Hz)')
ylabel('ab3A firing rate (Hz)')

% vary n, measure slopes and r^2
all_n = 1:21;
r2 = NaN(length(all_n),1);
slopes_lo = NaN(length(all_n),1);
slopes_hi = NaN(length(all_n),1);

if exist('.cache/variance_adapt_n.mat','file')
	load('.cache/variance_adapt_n.mat')
else
	for i = 1:length(all_n)
		textbar(i,length(all_n))
		q = p;
		q.n = all_n(i);

		Rhat = aNLN2(S,q);

		r2(i) = rsquare(VSdata.fA(1e4:end-1e4,example_trial), Rhat(1e4:end-1e4));

		Rhat = reshape(Rhat,1e4,length(Rhat)/1e4);
		Rhat(:,1) = []; Rhat(:,end) = [];

		ff = fit(vectorise(Rhat(1e3:5e3,:)),vectorise(R(1e3:5e3,:)),'poly1');
		slopes_hi(i) = ff.p1;

		ff = fit(vectorise(Rhat(6e3:end,:)),vectorise(R(6e3:end,:)),'poly1');
		slopes_lo(i) = ff.p1;
		
	end
	save('.cache/variance_adapt_n.mat','r2','slopes_lo','slopes_hi')
end

subplot(2,4,5); hold on
plot(all_n,slopes_lo(:)./slopes_hi(:),'k')
xlabel('n')
ylabel('gain_{low}/gain_{high}')
plot(all_n,0*all_n+1,'k--')

subplot(2,4,6); hold on
plot(all_n,r2,'k+-')
xlabel('n')
ylabel('r^2')

% now vary the adaptation parameters 
all_tau_z = unique(round(logspace(1,4,31)));
all_B = unique(logspace(-1,1,31));
r2 = NaN(length(all_tau_z),length(all_B));
slopes_hi = NaN(length(all_tau_z),length(all_B));
slopes_lo = NaN(length(all_tau_z),length(all_B));
if exist('.cache/variance_adapt_tau_B.mat','file')
	load('.cache/variance_adapt_tau_B.mat')
else
	for i = 1:length(all_tau_z)
		textbar(i,length(all_tau_z))
		for j = 1:length(all_B)
			q = p;
			q.tau_z = all_tau_z(i);
			q.B = all_B(j);

			Rhat = aNLN2(S,q);

			r2(i,j) = rsquare(VSdata.fA(1e4:end-1e4,example_trial), Rhat(1e4:end-1e4));

			Rhat = reshape(Rhat,1e4,length(Rhat)/1e4);
			Rhat(:,1) = []; Rhat(:,end) = [];

			ff = fit(vectorise(Rhat(1e3:5e3,:)),vectorise(R(1e3:5e3,:)),'poly1');
			slopes_hi(i,j) = ff.p1;

			ff = fit(vectorise(Rhat(6e3:end,:)),vectorise(R(6e3:end,:)),'poly1');
			slopes_lo(i,j) = ff.p1;

		end
	end
	save('.cache/variance_adapt_tau_B.mat','r2','slopes_lo','slopes_hi')
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

ax = subplot(2,4,7); hold on
imagesc(r2)
ylabel('\tau_{adaptation} (ms)')
xlabel('\beta')
h = colorbar;
caxis([0 1])
set(gca,'XTick',find(b_tick),'XTickLabel',B_labels(b_tick),'XTickLabelRotation',45)
set(gca,'YTick',find(tau_tick),'YTickLabel',tau_labels(tau_tick),'YTickLabelRotation',45)
set(gca,'XLim',[.5 size(r2,2)+.5],'YLim',[.5 size(r2,1)+.5])
title('r^2(data, adapting NLN model)')

ax2 = subplot(2,4,8); hold on
imagesc(abs(log(slopes_lo./slopes_hi)))
ylabel('\tau_{adaptation} (ms)')
xlabel('\beta')
colorbar
colormap(ax2,flipud(colormap))
caxis([0 10])
set(gca,'XTick',find(b_tick),'XTickLabel',B_labels(b_tick),'XTickLabelRotation',45)
set(gca,'YTick',find(tau_tick),'YTickLabel',tau_labels(tau_tick),'YTickLabelRotation',45)
set(gca,'XLim',[.5 size(r2,2)+.5],'YLim',[.5 size(r2,1)+.5])
title('abs(log(gain_{low}/gain_{high})')



prettyFig();

if being_published
	snapnow
	delete(gcf)
end

;;     ;;    ;;;    ;;;;;;;;  ;;;;    ;;;    ;;    ;;  ;;;;;;  ;;;;;;;; 
;;     ;;   ;; ;;   ;;     ;;  ;;    ;; ;;   ;;;   ;; ;;    ;; ;;       
;;     ;;  ;;   ;;  ;;     ;;  ;;   ;;   ;;  ;;;;  ;; ;;       ;;       
;;     ;; ;;     ;; ;;;;;;;;   ;;  ;;     ;; ;; ;; ;; ;;       ;;;;;;   
 ;;   ;;  ;;;;;;;;; ;;   ;;    ;;  ;;;;;;;;; ;;  ;;;; ;;       ;;       
  ;; ;;   ;;     ;; ;;    ;;   ;;  ;;     ;; ;;   ;;; ;;    ;; ;;       
   ;;;    ;;     ;; ;;     ;; ;;;; ;;     ;; ;;    ;;  ;;;;;;  ;;;;;;;; 

;;       ;;;;;;;; ;;;;;;;;  
;;       ;;       ;;     ;; 
;;       ;;       ;;     ;; 
;;       ;;;;;;   ;;;;;;;;  
;;       ;;       ;;        
;;       ;;       ;;        
;;;;;;;; ;;       ;;        


%%
% This data seems to strongly constrain the $n$, or the co-operativity parameter of the input Hill function. However, we know from other analysis that only part of the variance gain control occurs at the spiking machinery, and a equal portion occurs at transduction. Thus, if we are to interpret $n$ literally, I need to repeat this fit on the LFP data, as that is our proxy for receptor activity, not the firing rate. 

%%
% In the following figure, I vary the steepness of the input nonlinearity, and see how well it can account for the observed gain change in the LFP. All data is mean with s.e.m error bars. It looks like this data strongly constrains $n$ to be about 4. 

% remove baseline from all LFP 

PID = VSdata.PID;
LFP = VSdata.LFP;

global_start = 40e3; % 40 seconds
global_end = length(PID) - 5e3; 

% filter the LFP
for i = 1:width(LFP)
	LFP(:,i) = LFP(:,i) - fastFiltFilt(ones(1e4,1),1e4,LFP(:,i));
	LFP(:,i) = LFP(:,i)*10; % to get the units right, now in mV
end

% reshape the LFP signals
block_length = 1e4;
reshaped_LFP = LFP(global_start:end-1e4-1,1:width(PID));
reshaped_LFP = reshape(reshaped_LFP,block_length,width(reshaped_LFP)*length(reshaped_LFP)/block_length);

% also reshape the PID
reshaped_PID = PID(global_start:end-1e4-1,1:width(PID));
reshaped_PID = reshape(reshaped_PID,block_length,width(reshaped_PID)*length(reshaped_PID)/block_length);

% throw our NaNs globally. so we're throwing out epochs where the data is incomplete
rm_this = isnan(sum(reshaped_LFP));
reshaped_LFP(:,rm_this) = [];
reshaped_PID(:,rm_this) = [];

% filter to remove spikes
for i = 1:width(reshaped_LFP)
	reshaped_LFP(:,i) = filtfilt(ones(30,1),30,reshaped_LFP(:,i));
end

% project using filter
ft = 1e-3*(1:length(VSdata.K1)) - .1;
for i = 1:width(reshaped_PID)
	K1p(:,i) = convolve(1e-3*(1:length(reshaped_PID)),reshaped_PID(:,i),VSdata.K1,ft);
end




% calcualte gains while fitting a input nonlinearity -- trialwise
if exist('.cache/one_model_variance_NL_gain.mat','file') > 0
else
	all_n = 0:8;
	gain_lo = NaN(size(reshaped_PID,2),length(all_n));
	gain_lo2 = NaN(size(reshaped_PID,2),length(all_n));
	gain_hi = NaN(size(reshaped_PID,2),length(all_n));
	gain_hi2 = NaN(size(reshaped_PID,2),length(all_n));
	r2_lo = NaN(size(reshaped_PID,2),length(all_n));
	r2_hi = NaN(size(reshaped_PID,2),length(all_n));
	r2 = NaN(size(reshaped_PID,2),length(all_n));
	for i = 1:size(reshaped_PID,2)
		textbar(i,size(reshaped_PID,2))
		S = reshaped_PID(:,i);
		R = reshaped_LFP(:,i);
		k_D = mean(reshaped_PID(1e3:end-1e3));
		for j = 1:length(all_n)
			n = all_n(j);
			if n > 0
				X = 1./(1 + (k_D./S).^n);
				X = convolve(1e-3*(1:length(X)),X,VSdata.K1,ft);
				[ff,gof] = fit(K1p(1e3:5e3,i),X(1e3:5e3),'poly1');
				r2_hi(i,j) = gof.rsquare;
				gain_hi(i,j) = ff.p1;
				[ff,gof] = fit(K1p(6e3:9.5e3,i),X(6e3:9.5e3),'poly1');
				r2_lo(i,j) = gof.rsquare;
				gain_lo(i,j) = ff.p1;
				r2(i,j) = rsquare(K1p(1e3:end-1e3,i),X(1e3:end-1e3));

				% also estimate gain using the std. dev. ratios
				gain_lo2(i,j) = nanstd(X(6e3:9e3))/nanstd(S(6e3:9e3));
				gain_hi2(i,j) = nanstd(X(1e3:5e3))/nanstd(S(1e3:5e3));

			else
				[ff,gof]  = fit(K1p(1e3:5e3,i),R(1e3:5e3),'poly1');
				r2_hi(i,j) = gof.rsquare;
				gain_hi(i,j) = ff.p1;
				[ff,gof]  = fit(K1p(6e3:9.5e3,i),R(6e3:9.5e3),'poly1');
				gain_lo(i,j) = ff.p1;
				r2_lo(i,j) = gof.rsquare;
				r2(i,j) = rsquare(K1p(1e3:end-1e3,i),R(1e3:end-1e3));
			end
		end
	end
	save('.cache/one_model_variance_NL_gain.mat','all_n','gain_lo','gain_hi','r2','r2_lo','r2_hi');
end

% calcualte gains while fitting a input nonlinearity -- trialwise
if exist('.cache/one_model_variance_aNL_gain.mat','file') > 0

else
	all_n = 0:8;
	gain_lo = NaN(size(reshaped_PID,2),length(all_n));
	gain_hi = NaN(size(reshaped_PID,2),length(all_n));
	r2_lo = NaN(size(reshaped_PID,2),length(all_n));
	r2_hi = NaN(size(reshaped_PID,2),length(all_n));
	r2 = NaN(size(reshaped_PID,2),length(all_n));
	for i = 1:size(reshaped_PID,2)
		textbar(i,size(reshaped_PID,2))
		S = reshaped_PID(:,i);
		R = reshaped_LFP(:,i);
		k_D_hi = mean(reshaped_PID(1e3:5e3));
		k_D_lo = mean(reshaped_PID(6e3:end));
		for j = 1:length(all_n)
			n = all_n(j);
			if n > 0
				X = NaN*S;
				X(1:5e3) = 1./(1 + (k_D_hi./S(1:5e3)).^n);
				X(5e3+1:end) = 1./(1 + (k_D_lo./S(5e3+1:end)).^n);
				X = convolve(1e-3*(1:length(X)),X,VSdata.K1,ft);
				[ff,gof] = fit(K1p(1e3:5e3,i),X(1e3:5e3),'poly1');
				r2_hi(i,j) = gof.rsquare;
				gain_hi(i,j) = ff.p1;
				[ff,gof] = fit(K1p(6e3:9.5e3,i),X(6e3:9.5e3),'poly1');
				r2_lo(i,j) = gof.rsquare;
				gain_lo(i,j) = ff.p1;
				r2(i,j) = rsquare(K1p(1e3:end-1e3,i),X(1e3:end-1e3));
			else
				[ff,gof]  = fit(K1p(1e3:5e3,i),R(1e3:5e3),'poly1');
				r2_hi(i,j) = gof.rsquare;
				gain_hi(i,j) = ff.p1;
				[ff,gof]  = fit(K1p(6e3:9.5e3,i),R(6e3:9.5e3),'poly1');
				gain_lo(i,j) = ff.p1;
				r2_lo(i,j) = gof.rsquare;
				r2(i,j) = rsquare(K1p(1e3:end-1e3,i),R(1e3:end-1e3));
			end
		end
	end
	save('.cache/one_model_variance_aNL_gain.mat','all_n','gain_lo','gain_hi','r2','r2_lo','r2_hi','gain_hi2','gain_lo2');
end

% make a plot showing the range of gain change in the LFP, and which n accounts for this in the NL model
figure('outerposition',[0 0 1e3 500],'PaperUnits','points','PaperSize',[1000 500]); hold on


% first do the unadapting model
subplot(1,2,1); hold on
load('.cache/one_model_variance_NL_gain.mat')

temp = gain_lo(:,1)./gain_hi(:,1);
temp = sort(temp(r2_lo(:,1) > .5 & r2_hi(:,1) > .5));
m = mean(temp);
a = mean(temp) - std(temp)/sqrt(length(temp));
z = mean(temp) + std(temp)/sqrt(length(temp));

plot([0 8],[a a],'k--')
plot([0 8],[z z],'k--')
plot([0 8],[m m],'k','LineWidth',3)
set(gca,'YLim',[0.5 2])

y = NaN*(1:8);
ye = NaN*(1:8);
for i = 1:8
	temp = gain_lo(:,i+1)./gain_hi(:,i+1);
	temp = (temp(r2_lo(:,i+1) > .8 & r2_hi(:,i+1) > .8));
	y(i) = mean(temp);
	ye(i) = sem(temp);
end

errorbar(1:8,y,ye,'Color','r')

xlabel('n')
ylabel('gain_{low}/gain_{hi}')

title('NL model')


% now the model that accounts for the change in the mean
subplot(1,2,2); hold on
load('.cache/one_model_variance_aNL_gain.mat')

temp = gain_lo(:,1)./gain_hi(:,1);
temp = sort(temp(r2_lo(:,1) > .8 & r2_hi(:,1) > .8));
m = mean(temp);
a = mean(temp) - std(temp)/sqrt(length(temp));
z = mean(temp) + std(temp)/sqrt(length(temp));

plot([0 8],[a a],'k--')
plot([0 8],[z z],'k--')
plot([0 8],[m m],'k','LineWidth',3)
set(gca,'YLim',[0.5 2])

y = NaN*(1:8);
ye = NaN*(1:8);
for i = 1:8
	temp = gain_lo(:,i+1)./gain_hi(:,i+1);
	temp = (temp(r2_lo(:,i+1) > .5 & r2_hi(:,i+1) > .5));
	y(i) = mean(temp);
	ye(i) = sem(temp);
end

errorbar(1:8,y,ye,'Color','r')

xlabel('n')
ylabel('gain_{low}/gain_{hi}')

title('Accounting for change in \mu_{stim}')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end


;;;;;;;; ;;     ;; ;;       ;;          ;;       ;;;;;;;; ;;;;;;;;  
;;       ;;     ;; ;;       ;;          ;;       ;;       ;;     ;; 
;;       ;;     ;; ;;       ;;          ;;       ;;       ;;     ;; 
;;;;;;   ;;     ;; ;;       ;;          ;;       ;;;;;;   ;;;;;;;;  
;;       ;;     ;; ;;       ;;          ;;       ;;       ;;        
;;       ;;     ;; ;;       ;;          ;;       ;;       ;;        
;;        ;;;;;;;  ;;;;;;;; ;;;;;;;;    ;;;;;;;; ;;       ;;        

;;     ;;  ;;;;;;;  ;;;;;;;;  ;;;;;;;; ;;       
;;;   ;;; ;;     ;; ;;     ;; ;;       ;;       
;;;; ;;;; ;;     ;; ;;     ;; ;;       ;;       
;; ;;; ;; ;;     ;; ;;     ;; ;;;;;;   ;;       
;;     ;; ;;     ;; ;;     ;; ;;       ;;       
;;     ;; ;;     ;; ;;     ;; ;;       ;;       
;;     ;;  ;;;;;;;  ;;;;;;;;  ;;;;;;;; ;;;;;;;; 


%%
% In this section, I attempt to fit a kinetic model of the LFP (with binding/unbinding explicitly built in) to see what the $n$ parameter is here that accounts for the variance gain control we see in the data. 

clear data
c = 1;

for i = 1:10:size(reshaped_PID,2)
	data(c).stimulus = reshaped_PID(:,i);
	data(c).response = reshaped_LFP(:,i);
	data(c).response(1:1e3) = NaN;
	c = c+ 1;
end


%%
% First, I fit an adapting model that explicitly has binding kinetics to the data. The best-fit parameters of this model are:

clear p
p.       n = 5;
p.   k_min = 0.0112;
p. R_scale = -12.2080;
p.R_offset = -0.4607;
p.       B = 2.2062e+09;
p.      k2 = 3.0297;
p.     K_n = 0.1406;
p.   K_tau = 63.1562;

disp(p)


% calcualte gains and r2 
if exist('.cache/one_model_variance_LFPmodelv5F_gain.mat','file') > 0

else
	all_n = 0:8;
	gain_lo = NaN(size(reshaped_PID,2),length(all_n));
	gain_hi = NaN(size(reshaped_PID,2),length(all_n));
	r2_lo = NaN(size(reshaped_PID,2),length(all_n));
	r2_hi = NaN(size(reshaped_PID,2),length(all_n));
	r2 = NaN(size(reshaped_PID,2),length(all_n));
	for i = 1:size(reshaped_PID,2)
		textbar(i,size(reshaped_PID,2))
		S = reshaped_PID(:,i);
		R = reshaped_LFP(:,i);
		for j = 1:length(all_n)
			n = all_n(j);
			p.n = n;
			if n > 0
				X = vectorise(LFPmodelv5F(S,p));
				if ~any(isnan(X))
					[ff,gof] = fit(K1p(1e3:5e3,i),X(1e3:5e3),'poly1');
					r2_hi(i,j) = gof.rsquare;
					gain_hi(i,j) = ff.p1;
					[ff,gof] = fit(K1p(6e3:9.5e3,i),X(6e3:9.5e3),'poly1');
					r2_lo(i,j) = gof.rsquare;
					gain_lo(i,j) = ff.p1;
					r2(i,j) = rsquare(R(1e3:end-1e3),X(1e3:end-1e3));
				end
			else
			end
		end
	end
	save('.cache/one_model_variance_LFPmodelv5F_gain.mat','all_n','gain_lo','gain_hi','r2','r2_lo','r2_hi');
end


%%
% The $n$ parameter is 5, close to what the earlier, simpler analysis reported. What if I force the $n$ to be 1, and re-fit the model? The following figure shows the quality of predictions for the originally fit model, and a new model where I allow all paramters to vary, excpet $n$ which is fixed at 1. Note that the model with $n=5$ performs much better than the model where $n$ is nailed at 1. 

load('.cache/one_model_variance_LFPmodelv5F_gain.mat')
r2_n5 = r2(:,6);

clear p
p.       n = 1;
p.   k_min = 0.0737;
p. R_scale = -17.8377;
p.R_offset = -0.5876;
p.       B = 2.1392e+09;
p.      k2 = 674.5297;
p.     K_n = 1;
p.   K_tau = 85.2766;

% compute r2
r2_n1 = NaN*r2_n5;
for i = 1:length(r2_n1)
	S = reshaped_PID(:,i);
	R = reshaped_LFP(:,i);
	X = vectorise(LFPmodelv5F(S,p));
	r2_n1(i) = rsquare(X,R);
end


figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(r2_n1,r2_n5,'k+')
xlabel('r^2_{n=1}')
ylabel('r^2_{n=5}')
plot([0 1],[0 1],'k--')
set(gca,'XLim',[0 1],'YLim',[0 1])
prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%% Version Info
%
pFooter;


