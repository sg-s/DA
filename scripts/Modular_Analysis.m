% 
% 
% created by Srinivas Gorur-Shandilya at 4:36 , 21 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


pHeader;

%% Modular Analysis
% In this document, we take seriously the hypothesis that the LFP depends only on the stimulus, and the firing rate depends on the LFP, and analyze the data in these modules. Furthermore, we assume that the time derivative of the LFP is proportional to the transduction current. Thus, the (time-derivative) of the LFP should be predicted from the stimulus, and the firing rate should be predicted from the time-derivative of the LFP. 


%% Gain control by the mean stimulus
% In this section, we analyse the LFP to firing rate transformation in the mean shifted gaussian data. In the following figure, we extract filters from the PID, the derivative of the LFP and the firing rate, and compare them. 



[PID, LFP, fA, paradigm,~, ~, AllControlParadigms] = consolidateData('/local-data/DA-paper/LFP-MSG/september',1);


% sort the paradigms sensibly
sort_value = [];
for i = 1:length(AllControlParadigms)
	sort_value(i) = (mean(AllControlParadigms(i).Outputs(1,:)));
end
[~,idx] = sort(sort_value);

AllControlParadigms = AllControlParadigms(idx);
paradigm_new = paradigm*NaN;
for i = 1:length(idx)
	paradigm_new(paradigm == idx(i)) = i;
end
paradigm = paradigm_new;

% remove baseline from all PIDs
for i = 1:width(PID)
	PID(:,i) = PID(:,i) - mean(PID(1:5e3,i));
end

% throw out trials where we didn't record the LFP, for whatever reason
not_LFP = find((max(abs(LFP))) < 0.1);
LFP(:,not_LFP) = NaN;

% throw our bad traces
bad_trials = (sum(fA) == 0 | isnan(sum(fA)) |  isnan(sum(LFP)));
LFP(:,bad_trials) = [];
PID(:,bad_trials) = [];
fA(:,bad_trials) = [];
paradigm(bad_trials) = [];

% differentiate the LFP
dLFP = LFP;
for i = 1:length(paradigm)
	% first high pass them to remove spikes
	dLFP(:,i) = bandPass(LFP(:,i),Inf,10);
	dLFP(:,i) = 1e4*filtfilt(ones(10,1),10,[0; diff(dLFP(:,i))]);
end

% extract filters and find gain
a = 35e3; z = 55e3;
[K1,K1p,K1_gain] = extractFilters(PID,dLFP,'use_cache',true,'a',a,'z',z);
[K2,K2p,K2_gain] = extractFilters(dLFP,fA,'use_cache',true,'a',a,'z',z);
[K3,K3p,K3_gain] = extractFilters(PID,fA,'use_cache',true,'a',a,'z',z);
ft = 1e-3*(1:length(K1)) - .1;

% remember this for later
weber_data.K1 = K1;
weber_data.K2 = K2;
weber_data.K3 = K3;
weber_data.paradigm = paradigm;

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1), hold on
plot(ft,nanmean(K1(:,paradigm==1),2));
title('PID \rightarrow dLFP')
xlabel('Filtertime (s)')
ylabel('K1')

subplot(1,3,2), hold on
plot(ft,nanmean(K2(:,paradigm==1),2));
title('dLFP \rightarrow Firing','interpreter','tex')
xlabel('Filtertime (s)')
ylabel('K2')

subplot(1,3,3), hold on
plot(ft,nanmean(K3(:,paradigm==1),2));
temp1 = nanmean(K1(:,paradigm==1),2);
temp2 = nanmean(K2(:,paradigm==1),2);
plot(ft,convolve(ft,temp1,temp2,ft));
title('PID \rightarrow Firing','interpreter','tex')
legend('K3','K1 \otimes K2')

xlabel('Filtertime (s)')


prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end


%%
% Now, we compare the projections using these filters with the actual data. 


ss = 10;
figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1), hold on
x = nanmean(K1p(a:ss:z,paradigm==1),2);
y = nanmean(dLFP(a:ss:z,paradigm==1),2);
plot(x,y,'.');
title('PID \rightarrow dLFP')
xlabel('K1 \otimes s(t)')
ylabel('dLFP (mV/s)')
legend(['r^2 = ' oval(rsquare(x,y))],'Location','southeast')

subplot(1,3,2), hold on
x = nanmean(K2p(a:ss:z,paradigm==1),2);
y = nanmean(fA(a:ss:z,paradigm==1),2);
plot(x,y,'.');
title('dLFP \rightarrow Firing')
xlabel('K2 \otimes dLFP(t)')
ylabel('Firing Rate (Hz)')
legend(['r^2 = ' oval(rsquare(x,y))],'Location','southeast')

subplot(1,3,3), hold on
x = nanmean(K3p(a:ss:z,paradigm==1),2);
y = nanmean(fA(a:ss:z,paradigm==1),2);
plot(x,y,'.');
title('PID \rightarrow Firing')
xlabel('K3 \otimes s(t)')
ylabel('Firing Rate (Hz)')
legend(['r^2 = ' oval(rsquare(x,y))],'Location','southeast')

prettyFig()
if being_published	
	snapnow	
	delete(gcf)
end

%% 
% Now, we plot how the input-output curves of the dLFP and the firing rate change as a function of the mean stimulus. 

figure('outerposition',[0 0 1200 700],'PaperUnits','points','PaperSize',[1200 700]); hold on
ss = 100;
c = parula(max(paradigm)+1);
mean_stim = nanmean(PID(a:z,:));
ms = [min(mean_stim) max(mean_stim)];

subplot(2,3,1), hold on
for i = 1:max(paradigm) % iterate over all paradigms 
	y = nanmean(dLFP(a:z,paradigm == i),2);
	x = nanmean(K1p(a:z,paradigm == i),2);
	plotPieceWiseLinear(x,y,'nbins',50,'Color',c(i,:));
end
xlabel('K1 \otimes s(t)')
ylabel('dLPF (mV/s)')

subplot(2,3,2), hold on
for i = 1:max(paradigm) % iterate over all paradigms 
	y = nanmean(fA(a:z,paradigm == i),2);
	x = nanmean(K2p(a:z,paradigm == i),2);
	plotPieceWiseLinear(x,y,'nbins',50,'Color',c(i,:));
end
xlabel('K2 \otimes dLFP(t)')
ylabel('Firing Rate (Hz)')

subplot(2,3,3), hold on
for i = 1:max(paradigm) % iterate over all paradigms 
	y = nanmean(fA(a:z,paradigm == i),2);
	x = nanmean(K3p(a:z,paradigm == i),2);
	x = x - nanmean(x);
	s = nanmean(PID(a:z,paradigm==i),2);
	x = x + nanmean(s);
	plotPieceWiseLinear(x,y,'nbins',50,'Color',c(i,:));
end
xlabel('K3 \otimes s(t)')
ylabel('Firing Rate (Hz)')

subplot(2,3,4), hold on
for i = 1:length(paradigm)
	plot(mean_stim(i),K1_gain(i),'+','Color',c(paradigm(i),:))
end
xlabel('Mean Stimulus (V)')
ylabel('dLFP Gain (mV/s/V)')
set(gca,'XScale','log','YScale','log')
ff = fit(mean_stim(:),K1_gain(:),'power1','Upper',[Inf -1],'Lower',[0 -1]);
plot(ms,ff(ms),'r')


subplot(2,3,5), hold on
for i = 1:length(paradigm)
	plot(mean_stim(i),K2_gain(i),'+','Color',c(paradigm(i),:))
end
xlabel('Mean Stimulus (V)')
ylabel('Firing Gain (Hz/mV/s)')
set(gca,'XScale','log','YScale','log','YLim',[.10 10])


subplot(2,3,6), hold on
for i = 1:length(paradigm)
	plot(mean_stim(i),K3_gain(i),'+','Color',c(paradigm(i),:))
end
xlabel('Mean Stimulus (V)')
ylabel('Firing Gain (Hz/V)')
set(gca,'XScale','log','YScale','log')
ff = fit(mean_stim(:),K3_gain(:),'power1','Upper',[Inf -1],'Lower',[0 -1]);
plot(ms,ff(ms),'r')
prettyFig('fs',18)

if being_published	
	snapnow	
	delete(gcf)
end

%%
% This clearly demonstrates that gain control to the mean happens at the transduction level, which is responsible for the Weber-Fechner Law at the firing machinery. We can also see that the gain of the firing machinery itself doesn't change with the mean stimulus. 

%% 
% In this section, we try to account for the observed Weber-like scaling of transduction currents using a simple divisive gain control term:
% 
% $$ g(t)=\frac{\alpha}{1+\beta K_{\mu}\otimes s(t)} $$
% 

%%
% The best-fit model parameters were: 


p.tau = 79.8562;
p.  n = 1.0002;
p.  B = 4.2212;
p.  A = 213.57;

disp(p)

% remember this for later
weber_data.p = p;

%%
% However, the filter was very poorly constrained, and many different filters with lengths from 30ms to 3s seemed to do as good a job. 

d.stimulus = [PID(a:z,:), K1p(a:z,:)];
G = divisiveGainModel(d.stimulus,p);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
plot(mean_stim,K1_gain,'k+')
plot(mean_stim,G,'ro')
xlabel('Mean Stimulus')
legend('Measured','Model Prediction')
ylabel('Gain (mV/s/V)')
set(gca,'XScale','log','YScale','log')

subplot(1,2,2), hold on
plot([10 200],[10 200],'k--')
plot(G,K1_gain,'+')
xlabel('Model Prediction')
ylabel('Measured Gain (mV/s/V)')
set(gca,'XScale','log','YScale','log')

prettyFig('fs',18)

if being_published	
	snapnow	
	delete(gcf)
end

%% Contrast Sensitivity 
% In this section we analyse the data where we switch between two contrasts, while keeping the mean stimulus the same. First, we show the filters extracted at each stage. 


path_name = '/local-data/DA-paper/switching/variance/v2/';
[PID, LFP, fA, ~, orn] = consolidateData(path_name,1);

global_start = 40e3; % 40 seconds
global_end = length(PID) - 5e3; 
% bandpass to remove spikes and slow fluctuations
for i = 1:width(LFP)
	LFP(:,i) = bandPass(LFP(:,i),Inf,10); % remove spikes
	LFP(:,i) = 1e4*filtfilt(ones(10,1),10,[0; diff(LFP(:,i))]);
end

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

% extract filters and find gain
a = 1; z = 10e3;
K1 = extractFilters(reshaped_PID,reshaped_LFP,'use_cache',true,'a',a,'z',z);
K2 = extractFilters(reshaped_LFP,reshaped_fA,'use_cache',true,'a',a,'z',z);
K3 = extractFilters(reshaped_PID,reshaped_fA,'use_cache',true,'a',a,'z',z);
ft = 1e-3*(1:length(K1)) - .1;

% average filters and project stimulus
K1 = nanmean(K1,2);
K2 = nanmean(K2,2);
K3 = nanmean(K3,2);

K1p = NaN*reshaped_fA;
K2p = NaN*reshaped_fA;
K3p = NaN*reshaped_fA;
for i = 1:width(reshaped_fA)
	K1p(:,i) = convolve(1e-3*(1:length(reshaped_PID)),reshaped_PID(:,i),K1,ft);
	K2p(:,i) = convolve(1e-3*(1:length(reshaped_PID)),K1p(:,i),K2,ft);
	K3p(:,i) = convolve(1e-3*(1:length(reshaped_PID)),reshaped_PID(:,i),K3,ft);
end


figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1), hold on
plot(ft,K1);
title('PID \rightarrow dLFP')
xlabel('Filtertime (s)')
ylabel('K1')

subplot(1,3,2), hold on
plot(ft,K2);
title('dLFP \rightarrow Firing','interpreter','tex')
xlabel('Filtertime (s)')
ylabel('K2')

subplot(1,3,3), hold on
plot(ft,K3);
temp1 = K1;
temp2 = K2;
temp = convolve(ft,temp1,temp2,ft);
plot(ft,temp);
title('PID \rightarrow Firing','interpreter','tex')
legend('K3','K1 \otimes K2')
xlabel('Filtertime (s)')

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end


%%
% Now, we compare the projections using these filters with the actual data. 


ss = 10;
figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1), hold on
x = K1p(1e3:5e3,:);
y = reshaped_LFP(1e3:5e3,:);
rm_this = (isnan(sum(y)) | isnan(sum(x)));
x(:,rm_this) = []; y(:,rm_this) = []; 
h = plotPieceWiseLinear(x,y,'nbins',50,'Color',[1 0 0]);
delete(h.line(2:3))
delete(h.shade)
x = K1p(6e3:end-100,:);
y = reshaped_LFP(6e3:end-100,:);
rm_this = (isnan(sum(y)) | isnan(sum(x)));
x(:,rm_this) = []; y(:,rm_this) = []; 
h = plotPieceWiseLinear(x,y,'nbins',50,'Color',[0 0 1]);
delete(h.line(2:3))
delete(h.shade)
xlabel('K1 \otimes s(t)')
ylabel('dLPF (mV/s)')

subplot(1,3,2), hold on
x = K2p(1e3:5e3,:);
y = reshaped_fA(1e3:5e3,:);
rm_this = (isnan(sum(y)) | isnan(sum(x)));
x(:,rm_this) = []; y(:,rm_this) = []; 
h = plotPieceWiseLinear(x,y,'nbins',50,'Color',[1 0 0]);
delete(h.line(2:3))
delete(h.shade)
x = K2p(6e3:end-200,:);
y = reshaped_fA(6e3:end-200,:);
rm_this = (isnan(sum(y)) | isnan(sum(x)));
x(:,rm_this) = []; y(:,rm_this) = []; 
h = plotPieceWiseLinear(x,y,'nbins',50,'Color',[0 0 1]);
delete(h.line(2:3))
delete(h.shade)
xlabel('K2 \otimes dLFP_{pred}(t)')
ylabel('Firing Rate (Hz)')

subplot(1,3,3), hold on
x = K3p(1e3:5e3,:);
y = reshaped_fA(1e3:5e3,:);
rm_this = (isnan(sum(y)) | isnan(sum(x)));
x(:,rm_this) = []; y(:,rm_this) = []; 
h = plotPieceWiseLinear(x,y,'nbins',50,'Color',[1 0 0]);
delete(h.line(2:3))
delete(h.shade)
x = K3p(6e3:end-100,:);
y = reshaped_fA(6e3:end-100,:);
rm_this = (isnan(sum(y)) | isnan(sum(x)));
x(:,rm_this) = []; y(:,rm_this) = []; 
h = plotPieceWiseLinear(x,y,'nbins',50,'Color',[0 0 1]);
delete(h.line(2:3))
delete(h.shade)
xlabel('K3 \otimes s(t)')
ylabel('Firing Rate (Hz)')
prettyFig;

if being_published	
	snapnow	
	delete(gcf)
end

%% Naturalistic Stimuli: LFP
% In this section, we attempt to explain the LFP responses of neurons to naturalistic stimuli using what we learnt from the Weber experiment. 

[PID, LFP, fA, paradigm, orn, ~, AllControlParadigms] = consolidateData('/local-data/DA-paper/natural-flickering/with-lfp/ab3/',1);


% remove baseline from all PIDs
for i = 1:width(PID)
	PID(:,i) = PID(:,i) - mean(PID(1:5e3,i));
end

% throw out trials where we didn't record the LFP, for whatever reason
not_LFP = find((max(abs(LFP))) < 0.1);
LFP(:,not_LFP) = NaN;

% throw our bad traces
bad_trials = (sum(fA) == 0 | isnan(sum(fA)) |  isnan(sum(LFP)));
LFP(:,bad_trials) = [];
PID(:,bad_trials) = [];
fA(:,bad_trials) = [];
paradigm(bad_trials) = [];
orn(bad_trials) = [];

% differentiate the LFP
dLFP = LFP;
for i = 1:length(paradigm)
	% first high pass them to remove spikes
	dLFP(:,i) = 10*bandPass(LFP(:,i),Inf,30);
	dLFP(:,i) = 1e3*filtfilt(ones(30,1),30,[0; diff(dLFP(:,i))]);
end

% average across neurons
temp_PID = zeros(length(dLFP),max(orn));
temp_LFP = zeros(length(dLFP),max(orn));
temp_fA = zeros(length(dLFP),max(orn));
for i = 1:max(orn)
	temp_PID(:,i) = nanmean(PID(:,orn==i),2);
	temp_LFP(:,i) = nanmean(dLFP(:,orn==i),2);
	temp_fA(:,i) = nanmean(fA(:,orn==i),2);
end
PID = temp_PID; clear temp_PID
LFP = temp_LFP; clear temp_LFP dLFP
fA = temp_fA; clear temp_fA 
clear fly orn paradigm

% extract filters and find gain
a = 10e3; z = 60e3;
[K1,K1p,K1_gain] = extractFilters(PID,LFP,'use_cache',true,'a',a,'z',z);
[K2,K2p,K2_gain] = extractFilters(LFP,fA,'use_cache',true,'a',a,'z',z);
[K3,K3p,K3_gain] = extractFilters(PID,fA,'use_cache',true,'a',a,'z',z);
ft = 1e-3*(1:length(K1)) - .1;

%%
% In the following figure, we compare the LFP filter extracted from the natural stimuli with the LFP filter extracted from the Weber experiment, and also compare their linear projections, before and after gain-correction from the Weber-Fechner experiment.  

example_orn = 4;
figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1), hold on
plot(ft,K1(:,example_orn),'k')
plot(ft,nanmean(weber_data.K1(:,weber_data.paradigm == 1),2),'r')
legend('Natural Stimuli','Weber')
title('PID \rightarrow dLFP')
xlabel('Filtertime (s)')
ylabel('K1')

subplot(1,3,2), hold on
x = K1p(a:z,example_orn);
y = LFP(a:z,example_orn);
plot(x,y,'.')
legend(['r^2 = ' oval(rsquare(x,y))],'Location','southeast')
xlabel('K1 \otimes s(t)')
ylabel('dLFP (mV/s)')
title('No gain correction')

subplot(1,3,3), hold on
S = [K1p(:,example_orn), PID(:,example_orn)];
R = adaptiveGainModel(S,weber_data.p);
plot(R(a:z),LFP(a:z,example_orn),'.')
legend(['r^2 = ' oval(rsquare(R(a:z),LFP(a:z,example_orn)))],'Location','southeast')
xlabel('g(t)(K1 \otimes s(t))')
title('Corrected by Weber-calculated Gain')

prettyFig;

if being_published	
	snapnow	
	delete(gcf)
end



%% Version Info
%
pFooter;


