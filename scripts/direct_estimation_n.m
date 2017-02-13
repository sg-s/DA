
pHeader;


%% Estimating n directly from the data
% In this document I collect all the data we have, and plot some statistics of the data to get a sense of it. The ultimate goal is to try to directly infer $n$ from the data, without a model. 

figure('outerposition',[0 0 1001 901],'PaperUnits','points','PaperSize',[1001 901]); hold on
clear ax
for i = 1:4
	ax(i) = subplot(2,2,i); hold on
end

clear l


clear MSGdata
MSGdata = consolidateData2(getPath(dataManager,'93ba5d68174e3df9f462a1fc48c581da'));
MSGdata = cleanMSGdata(MSGdata);


Rg = NaN*MSGdata.paradigm;
Sg = NaN*MSGdata.paradigm;
for i = 1:length(MSGdata.paradigm)
	Rg(i) = std(MSGdata.fA(35e3:55e3,i))/mean(MSGdata.fA(35e3:55e3,i));
	Sg(i) = std(MSGdata.PID(35e3:55e3,i))/mean(MSGdata.PID(35e3:55e3,i));
end
mean_stim = nanmean(MSGdata.PID(35e3:55e3,:));
mean_resp = nanmean(MSGdata.fA(35e3:55e3,:));
std_resp = nanstd(MSGdata.PID(35e3:55e3,:));
l(1) = plot(ax(1),mean_stim,Rg./Sg,'k+');
plot(ax(2),mean_resp,Rg./Sg,'k+')
plot(ax(3),std_resp,Rg./Sg,'k+')

% now plot the same thing vs. the auto-correlation time of the stimulus
S_tau = autoCorrelationTime(MSGdata.PID(35e3:55e3,:));
plot(ax(4),S_tau,Rg./Sg,'k+');


;;     ;;    ;;;    ;;;;;;;;  ;;;;    ;;;    ;;    ;;  ;;;;;;  ;;;;;;;; 
;;     ;;   ;; ;;   ;;     ;;  ;;    ;; ;;   ;;;   ;; ;;    ;; ;;       
;;     ;;  ;;   ;;  ;;     ;;  ;;   ;;   ;;  ;;;;  ;; ;;       ;;       
;;     ;; ;;     ;; ;;;;;;;;   ;;  ;;     ;; ;; ;; ;; ;;       ;;;;;;   
 ;;   ;;  ;;;;;;;;; ;;   ;;    ;;  ;;;;;;;;; ;;  ;;;; ;;       ;;       
  ;; ;;   ;;     ;; ;;    ;;   ;;  ;;     ;; ;;   ;;; ;;    ;; ;;       
   ;;;    ;;     ;; ;;     ;; ;;;; ;;     ;; ;;    ;;  ;;;;;;  ;;;;;;;; 


clear VSdata
[VSdata.PID, VSdata.LFP, VSdata.fA, VSdata.paradigm, VSdata.orn, VSdata.fly] = consolidateData(getPath(dataManager,'e30707e8e8ef6c0d832eee31eaa585aa'),1);
% remove baseline from stimulus
VSdata.PID = bsxfun(@minus, VSdata.PID, min(VSdata.PID));
VSdata.LFP = bsxfun(@minus, VSdata.LFP, mean(VSdata.LFP(1:4000,:)));
v2struct(VSdata)

global_start = 40e3; % 40 seconds
global_end = length(PID) - 5e3; 

% reshape the LFP signals
block_length = 1e4;
reshaped_LFP = LFP(global_start:end-1e4-1,1:width(PID));
reshaped_LFP = reshape(reshaped_LFP,block_length,width(reshaped_LFP)*length(reshaped_LFP)/block_length);

% also reshape the PID
reshaped_PID = PID(global_start:end-1e4-1,1:width(PID));
reshaped_PID = reshape(reshaped_PID,block_length,width(reshaped_PID)*length(reshaped_PID)/block_length);

% reshape the firing rate signals
reshaped_fA = fA(global_start:end-1e4-1,1:width(PID));
reshaped_fA = reshape(reshaped_fA,block_length,width(reshaped_fA)*length(reshaped_fA)/block_length);


% also reshape the orn ID
reshaped_orn = repmat(orn,length(global_start:length(PID)-1e4-1)/block_length,1);
reshaped_orn = reshaped_orn(:);

% throw our NaNs globally. so we're throwing out epochs where the data is incomplete
rm_this = isnan(sum(reshaped_LFP));
reshaped_LFP(:,rm_this) = [];
reshaped_PID(:,rm_this) = [];
reshaped_fA(:,rm_this) = [];
reshaped_orn(rm_this) = [];
% 

std_hi = std(reshaped_PID(1e3:5e3,:));
std_lo = std(reshaped_PID(6e3:end,:));
mean_hi = mean(reshaped_PID(1e3:5e3,:));
mean_lo = mean(reshaped_PID(6e3:end,:));

mean_hi_R = mean(reshaped_fA(1e3:5e3,:));
mean_lo_R = mean(reshaped_fA(6e3:end,:));



Rg_lo = NaN(size(reshaped_PID,2),1);
Sg_lo = NaN(size(reshaped_PID,2),1);
Rg_hi = NaN(size(reshaped_PID,2),1);
Sg_hi = NaN(size(reshaped_PID,2),1);
for i = 1:length(Rg_hi)
	Sg_hi(i) = std(reshaped_PID(1e3:5e3,i))/mean(reshaped_PID(1e3:5e3,i));
	Sg_lo(i) = std(reshaped_PID(6e3:10e3,i))/mean(reshaped_PID(6e3:10e3,i));
	Rg_hi(i) = std(reshaped_fA(1e3:5e3,i))/mean(reshaped_fA(1e3:5e3,i));
	Rg_lo(i) = std(reshaped_fA(6e3:10e3,i))/mean(reshaped_fA(6e3:10e3,i));

end
plot(ax(2),mean_lo_R,Rg_lo./Sg_lo,'b+')
plot(ax(2),mean_hi_R,Rg_hi./Sg_hi,'r+')

plot(ax(3),std_lo,Rg_lo./Sg_lo,'b+')
plot(ax(3),std_hi,Rg_hi./Sg_hi,'r+')


% also plot it on the first plot
l(2) = plot(ax(1),mean_lo,Rg_lo./Sg_lo,'b+');
l(3) = plot(ax(1),mean_hi,Rg_hi./Sg_hi,'r+');



S_tau_lo = autoCorrelationTime(reshaped_PID(6e3:9e3,:));
S_tau_hi = autoCorrelationTime(reshaped_PID(1e3:5e3,:));

plot(ax(4),S_tau_lo,Rg_lo./Sg_lo,'b+')
plot(ax(4),S_tau_hi,Rg_hi./Sg_hi,'r+')



 ;;;;;;  ;;;;;;;;     ;;;    ;;;;;;;;   ;;;;;;  ;;;;;;;; 
;;    ;; ;;     ;;   ;; ;;   ;;     ;; ;;    ;; ;;       
;;       ;;     ;;  ;;   ;;  ;;     ;; ;;       ;;       
 ;;;;;;  ;;;;;;;;  ;;     ;; ;;;;;;;;   ;;;;;;  ;;;;;;   
      ;; ;;        ;;;;;;;;; ;;   ;;         ;; ;;       
;;    ;; ;;        ;;     ;; ;;    ;;  ;;    ;; ;;       
 ;;;;;;  ;;        ;;     ;; ;;     ;;  ;;;;;;  ;;;;;;;; 


% also add the sparse nat stim
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


Rg = NaN*SNSdata.orn;
Sg = NaN*SNSdata.orn;
for i = 1:length(SNSdata.orn)
	Rg(i) = std(SNSdata.fA(35e3:55e3,i))/mean(SNSdata.fA(35e3:55e3,i));
	Sg(i) = std(SNSdata.PID(35e3:55e3,i))/mean(SNSdata.PID(35e3:55e3,i));
end
mean_stim = nanmean(SNSdata.PID(35e3:55e3,:));
mean_resp = nanmean(SNSdata.fA(35e3:55e3,:));
std_resp = nanstd(SNSdata.PID(35e3:55e3,:));
l(4) = plot(ax(1),mean_stim,Rg./Sg,'go');
plot(ax(2),mean_resp,Rg./Sg,'go')
plot(ax(3),std_resp,Rg./Sg,'go')
plot(ax(4),autoCorrelationTime(SNSdata.PID),Rg./Sg,'go')



;;;;;;;;  ;;;;;;;; ;;    ;;  ;;;;;;  ;;;;;;;; 
;;     ;; ;;       ;;;   ;; ;;    ;; ;;       
;;     ;; ;;       ;;;;  ;; ;;       ;;       
;;     ;; ;;;;;;   ;; ;; ;;  ;;;;;;  ;;;;;;   
;;     ;; ;;       ;;  ;;;;       ;; ;;       
;;     ;; ;;       ;;   ;;; ;;    ;; ;;       
;;;;;;;;  ;;;;;;;; ;;    ;;  ;;;;;;  ;;;;;;;; 


% dense nat. stim
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

Rg = NaN*DNSdata.orn;
Sg = NaN*DNSdata.orn;
for i = 1:length(DNSdata.orn)
	Rg(i) = std(DNSdata.fA(35e3:55e3,i))/mean(DNSdata.fA(35e3:55e3,i));
	Sg(i) = std(DNSdata.PID(35e3:55e3,i))/mean(DNSdata.PID(35e3:55e3,i));
end
mean_stim = nanmean(DNSdata.PID(35e3:55e3,:));
mean_resp = nanmean(DNSdata.fA(35e3:55e3,:));
std_resp = nanstd(DNSdata.PID(35e3:55e3,:));
l(5) = plot(ax(1),mean_stim,Rg./Sg,'gd');
plot(ax(2),mean_resp,Rg./Sg,'gd')
plot(ax(3),std_resp,Rg./Sg,'gd')
plot(ax(4),autoCorrelationTime(DNSdata.PID),Rg./Sg,'gd')



;;       ;;     ;; ;;;;;;;; 
;;       ;;     ;; ;;       
;;       ;;     ;; ;;       
;;       ;;     ;; ;;;;;;   
;;        ;;   ;;  ;;       
;;         ;; ;;   ;;       
;;;;;;;;    ;;;    ;;       

p = getPath(dataManager,'c8dc5353a75ce5fcdcfa13139c716bd8');
[PID, LFP, fA, paradigm, orn, AllControlParadigms, paradigm_hashes,sequence] = consolidateData(p,1);

% remove baseline
PID = PID - min(min(PID));

Rg = NaN*paradigm;
Sg = NaN*paradigm;
for i = 1:length(paradigm)
	Rg(i) = std(fA(35e3:55e3,i))/mean(fA(35e3:55e3,i));
	Sg(i) = std(PID(35e3:55e3,i))/mean(PID(35e3:55e3,i));
end
mean_stim = nanmean(PID(35e3:55e3,:));
mean_resp = nanmean(fA(35e3:55e3,:));
std_resp = nanstd(PID(35e3:55e3,:));
l(6) = plot(ax(1),mean_stim,Rg./Sg,'mx');
plot(ax(2),mean_resp,Rg./Sg,'mx')
plot(ax(3),std_resp,Rg./Sg,'mx')

plot(ax(4),autoCorrelationTime(PID),Rg./Sg,'mx')



set(ax(1),'YLim',[0 8],'XScale','linear')
xlabel(ax(1),'\mu_{S} (V)')
ylabel(ax(1),'(\sigma_{R}/\sigma_{S})*(\mu_{S}/\mu_{R})')

xlabel(ax(2),'\mu_{R} (Hz)')
ylabel(ax(2),'(\sigma_{R}/\sigma_{S})*(\mu_{S}/\mu_{R})')
set(ax(2),'XLim',[0 60])

xlabel(ax(3),'\sigma_S (V)')
ylabel(ax(3),'(\sigma_{R}/\sigma_{S})*(\mu_{S}/\mu_{R})')

xlabel(ax(4),'\tau_{S} (ms)')
ylabel(ax(4),'(\sigma_{R}/\sigma_{S})*(\mu_{S}/\mu_{R})')

legend(l,{'Changing \mu','High \sigma','Low \sigma','Sparse nat.','Dense nat.','Exp. Gaussian'})

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end




%% Version Info
%
pFooter;

