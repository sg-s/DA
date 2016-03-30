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
	dLFP(:,i) = -1e4*filtfilt(ones(10,1),10,[0; diff(dLFP(:,i))]);
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

iv = sum(nanmean(K1(:,paradigm==1),2))./sum(abs(nanmean(K1(:,paradigm==1),2)));
tx = ['$\frac{\intop K}{\intop\left|K\right|}$=' oval(iv)];
text(.3,10e-3,tx,'interpreter','latex')

subplot(1,3,2), hold on
plot(ft,nanmean(K2(:,paradigm==1),2));
title('dLFP \rightarrow Firing','interpreter','tex')
xlabel('Filtertime (s)')
ylabel('K2')

iv = sum(nanmean(K2(:,paradigm==1),2))./sum(abs(nanmean(K2(:,paradigm==1),2)));
tx = ['$\frac{\intop K}{\intop\left|K\right|}$=' oval(iv)];
text(.3,10e-3,tx,'interpreter','latex')


subplot(1,3,3), hold on
plot(ft,nanmean(K3(:,paradigm==1),2));
temp1 = nanmean(K1(:,paradigm==1),2);
temp2 = nanmean(K2(:,paradigm==1),2);
plot(ft,convolve(ft,temp1,temp2,ft));
title('PID \rightarrow Firing','interpreter','tex')
legend('K3','K1 \otimes K2')

iv = sum(nanmean(K3(:,paradigm==1),2))./sum(abs(nanmean(K3(:,paradigm==1),2)));
tx = ['$\frac{\intop K}{\intop\left|K\right|}$=' oval(iv)];
text(.3,10e-3,tx,'interpreter','latex')

xlabel('Filtertime (s)')

labelFigure
prettyFig

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
labelFigure
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

labelFigure
prettyFig('fs',18)

if being_published	
	snapnow	
	delete(gcf)
end

%%
% This clearly demonstrates that gain control to the mean happens at the transduction level, which is responsible for the Weber-Fechner Law at the firing machinery. We can also see that the gain of the firing machinery itself doesn't change with the mean stimulus. 

%% Accounting for mean-sensitive gain changes
% In this section, we try to account for the observed Weber-like scaling of transduction currents using a simple divisive gain control term:
% 
% $$ g(t)=\frac{\alpha}{1+\beta K_{\mu}\otimes s(t)} $$
% 

%%
% The best-fit model parameters were: 


p.tau = 79.8562;
p.  n = 1.0002;
p.  B = 100;
p.  A = 3e3;

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

labelFigure
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
dLFP = LFP;
for i = 1:width(LFP)
	LFP(:,i) = bandPass(LFP(:,i),Inf,30); % remove spikes
	dLFP(:,i) = -1e4*filtfilt(ones(10,1),10,[0; diff(LFP(:,i))]);
end

% reshape the LFP signals
block_length = 1e4;
reshaped_LFP = LFP(global_start:end-1e4-1,1:width(PID));
reshaped_LFP = reshape(reshaped_LFP,block_length,width(reshaped_LFP)*length(reshaped_LFP)/block_length);

reshaped_dLFP = dLFP(global_start:end-1e4-1,1:width(PID));
reshaped_dLFP = reshape(reshaped_dLFP,block_length,width(reshaped_dLFP)*length(reshaped_dLFP)/block_length);

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
rm_this = isnan(sum(reshaped_dLFP));
reshaped_LFP(:,rm_this) = [];
reshaped_dLFP(:,rm_this) = [];
reshaped_PID(:,rm_this) = [];
reshaped_fA(:,rm_this) = [];
reshaped_orn(rm_this) = [];

% extract filters and find gain
a = 1; z = 10e3;
K1 = extractFilters(reshaped_PID,reshaped_dLFP,'use_cache',true,'a',a,'z',z);
K2 = extractFilters(reshaped_dLFP,reshaped_fA,'use_cache',true,'a',a,'z',z);
K3 = extractFilters(reshaped_PID,reshaped_fA,'use_cache',true,'a',a,'z',z);
ft = 1e-3*(1:length(K1)) - .1;

% throw out traces where the dLFP traces are significantly far from 0
dmax = (nanmean(nanmean(reshaped_dLFP)) + 2*nanstd(nanmean(reshaped_dLFP)));
dmin = (nanmean(nanmean(reshaped_dLFP)) - 2*nanstd(nanmean(reshaped_dLFP)));
rm_this = nanmean(reshaped_dLFP) < dmin | nanmean(reshaped_dLFP) > dmax;

reshaped_LFP(:,rm_this) = [];
reshaped_dLFP(:,rm_this) = [];
reshaped_PID(:,rm_this) = [];
reshaped_fA(:,rm_this) = [];
reshaped_orn(rm_this) = [];

% remove mean from the LFP for each trial
for i = 1:width(reshaped_LFP)
	reshaped_LFP(:,i) =  reshaped_LFP(:,i) - nanmean(reshaped_LFP(:,i));
end

% average filters and project stimulus
K1 = nanmean(K1,2);
K2 = nanmean(K2,2); K2(K2<0) = 0;
K3 = nanmean(K3,2);

K1p = NaN*reshaped_fA;
K2p = NaN*reshaped_fA;
K3p = NaN*reshaped_fA;
K2K1p = NaN*reshaped_fA;
for i = 1:width(reshaped_fA)
	K1p(:,i) = convolve(1e-3*(1:length(reshaped_PID)),reshaped_PID(:,i),K1,ft);
	K2K1p(:,i) = convolve(1e-3*(1:length(reshaped_PID)),K1p(:,i),K2,ft);
	K2p(:,i) = convolve(1e-3*(1:length(reshaped_PID)),reshaped_dLFP(:,i),K2,ft);
	K3p(:,i) = convolve(1e-3*(1:length(reshaped_PID)),reshaped_PID(:,i),K3,ft);
end


figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1), hold on
plot(ft,K1);
title('PID \rightarrow dLFP')
xlabel('Filtertime (s)')
ylabel('K1')

iv = sum(K1)./sum(abs(K1));
tx = ['$\frac{\intop K}{\intop\left|K\right|}$=' oval(iv)];
text(.3,10e-3,tx,'interpreter','latex')

subplot(1,3,2), hold on
plot(ft,K2);
title('dLFP \rightarrow Firing','interpreter','tex')
xlabel('Filtertime (s)')
ylabel('K2')

iv = sum(K2)./sum(abs(K2));
tx = ['$\frac{\intop K}{\intop\left|K\right|}$=' oval(iv)];
text(.3,10e-3,tx,'interpreter','latex')

subplot(1,3,3), hold on
plot(ft,K3);
temp1 = K1;
temp2 = K2;
temp = convolve(ft,temp1,temp2,ft);
plot(ft,temp);
title('PID \rightarrow Firing','interpreter','tex')
legend('K3','K1 \otimes K2')
xlabel('Filtertime (s)')

iv = sum(K3)./sum(abs(K3));
tx = ['$\frac{\intop K}{\intop\left|K\right|}$=' oval(iv)];
text(.3,10e-3,tx,'interpreter','latex')

labelFigure
prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%%
% Why is the first filter not a perfect differentiator? We would expect it to be so since we expect the derivative of the LFP to have a mean of zero. 

%%
% Now, we compare the projections using these filters with the actual data. In the following figure, I plot responses (dLFP or firing rate) vs various projections using the filters we backed out. In each case, the red curve is from the high variance epoch, and the blue curve is from the low variance epoch. For each of the plots, I also quantify the degree of gain control as follows. First, I compute the gain in the low and high variance epochs. This ratio (low variance gain/high variance gain) is the first measure. Second, I also compute the standard deviation of the input in the low variance epoch, and divide by the standard deviation of the input during the high variance epoch. In each plot, the "input" is whatever is being convolved by the outermost filter. So, for example, in (c), the input is the predicted dLFP.


% compute the r2 for each trial
r2_K1p = NaN(width(reshaped_dLFP),1);
r2_K2p = NaN(width(reshaped_dLFP),1);
r2_K2K1p = NaN(width(reshaped_dLFP),1);
r2_K3p = NaN(width(reshaped_dLFP),1);
for i = 1:length(r2_K1p)
	r2_K1p(i) = rsquare(K1p(1e3:9e3,i),reshaped_dLFP(1e3:9e3,i));
	r2_K2p(i) = rsquare(K2p(1e3:9e3,i),reshaped_fA(1e3:9e3,i));
	r2_K2K1p(i) = rsquare(K2K1p(1e3:9e3,i),reshaped_fA(1e3:9e3,i));
	r2_K3p(i) = rsquare(K3p(1e3:9e3,i),reshaped_fA(1e3:9e3,i));
end

ss = 1;
min_r2 = .7; 

if ~exist('input_contrast_change_K1p','var')
	computeContrastsGains;
end

figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[800 800]); hold on
subplot(2,2,1), hold on
x = K1p(1e3:4e3,:);
y = reshaped_dLFP(1e3:4e3,:);
rm_this = (isnan(sum(y)) | isnan(sum(x)) | (r2_K1p < min_r2)');
x(:,rm_this) = []; y(:,rm_this) = []; 
x = x(:); y = y(:);
x = x(1:ss:end); y = y(1:ss:end);
plotPieceWiseLinear(x(:),y(:),'nbins',50,'Color',[1 0 0],'proportional_bins',true,'show_error',false);
x = K1p(6e3:9e3,:);
y = reshaped_dLFP(6e3:9e3,:);
rm_this = (isnan(sum(y)) | isnan(sum(x)) | (r2_K1p < min_r2)');
x(:,rm_this) = []; y(:,rm_this) = []; 
x = x(:); y = y(:);
x = x(1:ss:end); y = y(1:ss:end);
plotPieceWiseLinear(x(:),y(:),'nbins',50,'Color',[0 0 1],'proportional_bins',true,'show_error',false);
xlabel('K1 \otimes s(t)')
ylabel('dLPF (mV/s)')

tx = ['$\Delta gain$ = ' oval(nanmean(gain_change_K1p)) ,'$ \pm $' oval(sem(gain_change_K1p))];
text(-.2,20,tx,'interpreter','latex');
tx = ['$\Delta \sigma $= ' oval(nanmean(input_contrast_change_K1p)) ,'$ \pm $' oval(sem(input_contrast_change_K1p))];
text(-.2,15,tx,'interpreter','latex');

subplot(2,2,2), hold on
x = K2p(1e3:4e3,:);
y = reshaped_fA(1e3:4e3,:);
rm_this = (isnan(sum(y)) | isnan(sum(x)) | (r2_K2p < min_r2)');
x(:,rm_this) = []; y(:,rm_this) = []; 
x = x(:); y = y(:);
x = x(1:ss:end); y = y(1:ss:end);
plotPieceWiseLinear(x(:),y(:),'nbins',50,'Color',[1 0 0],'proportional_bins',true,'show_error',false);
x = K2p(6e3:9e3,:);
y = reshaped_fA(6e3:9e3,:);
rm_this = (isnan(sum(y)) | isnan(sum(x)) | (r2_K2p < min_r2)');
x(:,rm_this) = []; y(:,rm_this) = []; 
x = x(:); y = y(:);
x = x(1:ss:end); y = y(1:ss:end);
plotPieceWiseLinear(x(:),y(:),'nbins',50,'Color',[0 0 1],'proportional_bins',true,'show_error',false);
xlabel('K2 \otimes dLFP(t)')
ylabel('Firing Rate (Hz)')

tx = ['$\Delta gain $= ' oval(nanmean(gain_change_K2p)) ,'$ \pm $' oval(sem(gain_change_K2p))];
text(-.2,20,tx,'interpreter','latex');
tx = ['$\Delta \sigma $= ' oval(nanmean(input_contrast_change_K2p)) ,'$ \pm $' oval(sem(input_contrast_change_K2p))];
text(-.2,15,tx,'interpreter','latex');


subplot(2,2,3), hold on
x = K2K1p(1e3:4e3,:);
y = reshaped_fA(1e3:4e3,:);
rm_this = (isnan(sum(y)) | isnan(sum(x)) | (r2_K2K1p < min_r2)');
x(:,rm_this) = []; y(:,rm_this) = []; 
x = x(:); y = y(:);
x = x(1:ss:end); y = y(1:ss:end);
plotPieceWiseLinear(x(:),y(:),'nbins',50,'Color',[1 0 0],'proportional_bins',true,'show_error',false);
x = K2K1p(6e3:9e3,:);
y = reshaped_fA(6e3:9e3,:);
rm_this = (isnan(sum(y)) | isnan(sum(x)) | (r2_K2K1p < min_r2)');
x(:,rm_this) = []; y(:,rm_this) = []; 
x = x(:); y = y(:);
x = x(1:ss:end); y = y(1:ss:end);
plotPieceWiseLinear(x(:),y(:),'nbins',50,'Color',[0 0 1],'proportional_bins',true,'show_error',false);
xlabel('K2 \otimes K1 \otimes s(t)')
ylabel('Firing Rate (Hz)')

tx = ['$\Delta gain $= ' oval(nanmean(gain_change_K2K1p)) ,'$ \pm $' oval(sem(gain_change_K2K1p))];
text(-.2,60,tx,'interpreter','latex');
tx = ['$\Delta \sigma $= ' oval(nanmean(input_contrast_change_K2K1p)) ,'$ \pm $' oval(sem(input_contrast_change_K2K1p))];
text(-.2,52,tx,'interpreter','latex');

subplot(2,2,4), hold on
x = K3p(1e3:4e3,:);
y = reshaped_fA(1e3:4e3,:);
rm_this = (isnan(sum(y)) | isnan(sum(x)) | (r2_K3p < min_r2)');
x(:,rm_this) = []; y(:,rm_this) = []; 
x = x(:); y = y(:);
x = x(1:ss:end); y = y(1:ss:end);
plotPieceWiseLinear(x(:),y(:),'nbins',50,'Color',[1 0 0],'proportional_bins',true,'show_error',false);
x = K3p(6e3:9e3,:);
y = reshaped_fA(6e3:9e3,:);
rm_this = (isnan(sum(y)) | isnan(sum(x)) | (r2_K3p < min_r2)');
x(:,rm_this) = []; y(:,rm_this) = []; 
x = x(:); y = y(:);
x = x(1:ss:end); y = y(1:ss:end);
plotPieceWiseLinear(x(:),y(:),'nbins',50,'Color',[0 0 1],'proportional_bins',true,'show_error',false);
xlabel('K3 \otimes  s(t)')
ylabel('Firing Rate (Hz)')

tx = ['$\Delta gain$ = ' oval(nanmean(gain_change_K3p)) ,'$ \pm $' oval(sem(gain_change_K3p))];
text(.1,60,tx,'interpreter','latex');
tx = ['$\Delta \sigma $= ' oval(nanmean(input_contrast_change_K3p)) ,'$ \pm $' oval(sem(input_contrast_change_K3p))];
text(.1,52,tx,'interpreter','latex');

labelFigure
prettyFig;

if being_published	
	snapnow	
	delete(gcf)
end

%%
% This is weird. The contrast adaptation is strongest not from the LFP to the firing rate, as we would expect, but from the the predicted LFP to the firing rate. To look at this more closely, I plot predicted dLFP vs. the actual dLFP for the high and low variance epochs (a). I'm also plotting the distributions of the dLFPs during the two epochs (b) and the distributions of the predicted dLFPs (c).

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1), hold on
plot(nanmean(K1p(1e3:4e3,:)),nanmean(reshaped_dLFP(1e3:4e3,:)),'r+')
plot(nanmean(K1p(6e3:9e3,:)),nanmean(reshaped_dLFP(6e3:9e3,:)),'b+')
xlabel('<dLFP_{pred}>_{trials}')
ylabel('<dLFP (mV/s)>_{trials}')

subplot(1,3,2), hold on
plot([0 0],[0 .14],'k--')
x = vectorise(reshaped_dLFP(1e3:4e3,:));
[hy,hx] = histcounts(x,linspace(min(x),max(x),100));
hy = hy/sum(hy);
hx = hx(1:end-1) + mean(diff(hx));
plot(hx,hy,'r')
x = vectorise(reshaped_dLFP(6e3:9e3,:));
[hy,hx] = histcounts(x,linspace(min(x),max(x),100));
hy = hy/sum(hy);
hx = hx(1:end-1) + mean(diff(hx));
plot(hx,hy,'b')
xlabel('dLFP (mV/s)')
ylabel('Probability')
set(gca,'XLim',[-40 40])

subplot(1,3,3), hold on
plot([0 0],[0 .04],'k--')
x = vectorise(K1p(1e3:4e3,:));
[hy,hx] = histcounts(x,linspace(min(x),max(x),100));
hy = hy/sum(hy);
hx = hx(1:end-1) + mean(diff(hx));
plot(hx,hy,'r')
x = vectorise(K1p(6e3:9e3,:));
[hy,hx] = histcounts(x,linspace(min(x),max(x),100));
hy = hy/sum(hy);
hx = hx(1:end-1) + mean(diff(hx));
plot(hx,hy,'b')
xlabel('dLFP_{pred}')
ylabel('Probability')
set(gca,'XLim',[-.4 .5])

labelFigure
prettyFig;

if being_published	
	snapnow	
	delete(gcf)
end

%%
% There is something very weird going on. Specifically, the red distribution seems to be shifted to the left of zero, suggesting that the LFP is continuously increasing during the high variance epoch (This is because all reported dLFPs are actually -dLFPs). Does this make any sense? In the following figure, I plot the raw LFP traces, mean subtracted across each trial, to see if there is any epoch-specific trend. I also fit lines to the high (red) and low (blue) variance epochs, to show that there is an upward or downward trend. 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
temp = (nanmean(reshaped_LFP,2));
x = 1e3:4e3;
y = temp(x);
ff = fit(x(:),y(:),'poly1');
plot(temp,'k')
plot(x,ff(x),'r')
x = 6e3:9e3;
y = temp(x);
ff = fit(x(:),y(:),'poly1');
plot(x,ff(x),'b')
xlabel('Time (ms)')
ylabel('<LFP>_{trials} (\DeltamV)')
prettyFig;

if being_published	
	snapnow	
	delete(gcf)
end

%%
% OK, so now things make sense. We see an upward trend in the LFP during the high variance epoch, which explains the red distribution being pushed to the left, and a downward trend in the low variance epoch, which explains the blue distribution pushed towards the right. 

%%
% Why is the mean LFP higher in the high-contrast case than in the low-contrast case? The simplest explanation is because the stimulus has slightly lower mean when the variance is high (because of an error in the stimulus delivery). 

%%
% Does this invalidate our result? Not really, for two reasons. First, the unbiased stimulus to firing rate model clearly shows contrast adaptation, and there is no getting around that. Second, the decrease in stimulus mean during the high contrast epoch, which probably causes the deviation in the LFP, works against us, and decreases the size of the expected effect, not increases it. 

%% Inserting the Weber Model into the contrast data
% In this section, we insert the Weber-Fechner Model we measured using the last dataset into this one, because we know it exists, and plays a role here. 

K1p_weber = K1p;
K2K1p_weber = NaN*K1p_weber;
for i = 1:width(K1p)
	d.stimulus = [K1p(:,i),reshaped_PID(:,i)];
	K1p_weber(:,i) = adaptiveGainModel(d.stimulus,weber_data.p);
	% now convolve with K2
	K2K1p_weber(:,i) = convolve(1e-3*(1:length(reshaped_PID)),K1p_weber(:,i),K2,ft);
end

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1), hold on
x = K1p(1e3:4e3,:);
y = reshaped_dLFP(1e3:4e3,:);
rm_this = (isnan(sum(y)) | isnan(sum(x)));
x(:,rm_this) = []; y(:,rm_this) = []; 
h = plotPieceWiseLinear(x,y,'nbins',50,'Color',[1 0 0]);
delete(h.line(2:3))
delete(h.shade)
x = K1p(6e3:9e3,:);
y = reshaped_dLFP(6e3:9e3,:);
rm_this = (isnan(sum(y)) | isnan(sum(x)));
x(:,rm_this) = []; y(:,rm_this) = []; 
h = plotPieceWiseLinear(x,y,'nbins',50,'Color',[0 0 1]);
delete(h.line(2:3))
delete(h.shade)
xlabel('K1 \otimes s(t)')
ylabel('dLFP (mV/s)')

subplot(1,3,2), hold on
x = K1p_weber(1e3:4e3,:);
y = reshaped_dLFP(1e3:4e3,:);
rm_this = (isnan(sum(y)) | isnan(sum(x)));
x(:,rm_this) = []; y(:,rm_this) = []; 
h = plotPieceWiseLinear(x,y,'nbins',50,'Color',[1 0 0]);
delete(h.line(2:3))
delete(h.shade)
x = K1p_weber(6e3:9e3,:);
y = reshaped_dLFP(6e3:9e3,:);
rm_this = (isnan(sum(y)) | isnan(sum(x)));
x(:,rm_this) = []; y(:,rm_this) = []; 
h = plotPieceWiseLinear(x,y,'nbins',50,'Color',[0 0 1]);
delete(h.line(2:3))
delete(h.shade)
xlabel('g_{Weber}(t) * K1 \otimes s(t)')
ylabel('dLFP (mV/s)')

subplot(1,3,3), hold on
x = K2K1p_weber(1e3:4e3,:);
y = reshaped_fA(1e3:4e3,:);
rm_this = (isnan(sum(y)) | isnan(sum(x)));
x(:,rm_this) = []; y(:,rm_this) = []; 
h = plotPieceWiseLinear(x,y,'nbins',50,'Color',[1 0 0]);
delete(h.line(2:3))
delete(h.shade)
x = K2K1p_weber(6e3:9e3,:);
y = reshaped_fA(6e3:9e3,:);
rm_this = (isnan(sum(y)) | isnan(sum(x)));
x(:,rm_this) = []; y(:,rm_this) = []; 
h = plotPieceWiseLinear(x,y,'nbins',50,'Color',[0 0 1]);
delete(h.line(2:3))
delete(h.shade)
xlabel('K2 \otimes (g_{Weber}(t) * K1 \otimes s(t))')
ylabel('Firing Rate (Hz)')

labelFigure
prettyFig;

if being_published	
	snapnow	
	delete(gcf)
end


%% Accounting for the Contrast Gain Change
% Now, we want to account for this gain change. We use a logistic function to describe the shape of the input-output curve from the LFP to the firing rate:
% 
% $$ f(x)=\frac{A}{1+e^{-kx-x_{0}}} $$
%

%%
% where k controls the steepness of the curve. 

%%
% First, we fit this logistic function to the two input-output curves, allowing only the steepness parameter to vary between the blue and red curves.

x0 = -1.3245;

ft = fittype('logistic(x,A,k,x0)');
ff_lo = fit(data_lo.x(:),data_lo.y(:),ft,'StartPoint',[max(data_hi.y) 2 x0],'Lower',[max(data_hi.y) 1 x0],'Upper',[max(data_hi.y) 100 x0],'MaxIter',1e6);
ff_hi = fit(data_hi.x(:),data_hi.y(:),ft,'StartPoint',[max(data_hi.y) 2 0],'Lower',[max(data_hi.y) x0],'Upper',[max(data_hi.y) 100 x0],'MaxIter',1e6);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
plot(data_lo.x,data_lo.y,'b')
plot(data_hi.x,ff_lo(data_hi.x),'b--')
title(['k = ' oval(ff_lo.k)])
ylabel('ORN response (Hz)')
xlabel('Projected Stimulus')

subplot(1,2,2), hold on
plot(data_hi.x,data_hi.y,'r')
plot(data_hi.x,ff_hi(data_hi.x),'r--')
title(['k = ' oval(ff_hi.k)])
ylabel('ORN response (Hz)')
xlabel('Projected Stimulus')

labelFigure
prettyFig;

if being_published	
	snapnow	
	delete(gcf)
end

%%
% The steepness of the curve should depend on variance on the LFP signal. This way, the firing rate gain that depends on the contrast depends only on various parts of the LFP signal, making it modular. We assume
% 
% $$ k=\frac{k_{0}}{1+\beta_{\sigma}K_{\sigma}\otimes\left[l'\right]_{+}} $$
% 
% and find the best-fit parameters that account for the data. 
% 

clear d
these_trials = 50:10:200;
for i = length(these_trials):-1:1
	d(i).stimulus = reshaped_dLFP(:,these_trials(i));
	d(i).response = [ff_hi.k,ff_lo.k];
end

clear p
p. k0 = 17.9078;
p.  B = 0.0432;
p.  n = 1.7656;
p.tau = 97.6875;
p. s0 = 4.3347;
p. s1 = -30.4246;


% add some parameters
p.A = ff_hi.A;
p.x0 = ff_hi.x0;

XG = K2p;
for i = 1:width(K2p)
	XG(:,i) = contrastLogisticModel([(reshaped_dLFP(:,i)*p.s0 + p.s1) ,K2K1p(:,i)],p);
end

figure('outerposition',[0 0 900 800],'PaperUnits','points','PaperSize',[900 800]); hold on
subplot(2,2,1), hold on
plot(data_lo.x,data_lo.y,'b')
plot(data_hi.x,data_hi.y,'r')
xlabel('K \otimes s')
ylabel('Response (Hz)')

subplot(2,2,2), hold on
plotPieceWiseLinear(XG(1e3:5e3,:),reshaped_fA(1e3:5e3,:),'nbins',50,'Color','r');
plotPieceWiseLinear(XG(6e3:9e3,:),reshaped_fA(6e3:9e3,:),'nbins',50,'Color','b');
xlabel('$\hat{R}$','interpreter','latex')
ylabel('Response (Hz)')

subplot(2,2,3); hold on
plotPieceWiseLinear(K2K1p(1e3:5e3,:),XG(1e3:5e3,:),'nbins',50,'Color','r');
plotPieceWiseLinear(K2K1p(6e3:9e3,:),XG(6e3:9e3,:),'nbins',50,'Color','b');
ylabel('$\hat{R}$','interpreter','latex')
xlabel('K \otimes s')

subplot(2,2,4); hold on
r2_X = NaN(width(K2K1p),1);
r2_XG = r2_X;
for i = 1:width(K2K1p)
	fp = K2K1p([1e3:5e3 6e3:9e3],i); r = reshaped_fA([1e3:5e3 6e3:9e3],i);
	try
		r2_X(i) = rsquare(fp,r);
	catch
	end
	fp = XG([1e3:5e3 6e3:9e3],i);
	try
		r2_XG(i) = rsquare(fp,r);
	catch
	end
end
plot([0 1],[0 1],'k--')
plot(r2_X,r2_XG,'k+')
xlabel('r^2 Linear prediction')
ylabel('r^2 contrast-corrected')

labelFigure
prettyFig('fs',16)

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
	dLFP(:,i) = bandPass(LFP(:,i),Inf,30);
	dLFP(:,i) = -1e4*filtfilt(ones(10,1),10,[0; diff(dLFP(:,i))]);
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

example_orn = 2;
Weber_K = nanmean(weber_data.K1(:,weber_data.paradigm == 1),2);
figure('outerposition',[0 0 1500 800],'PaperUnits','points','PaperSize',[1500 800]); hold on
subplot(2,3,1), hold on
plot(ft,K1(:,example_orn),'k')
plot(ft,Weber_K,'r')
legend('Natural Stimuli','Weber')
title('PID \rightarrow dLFP')
xlabel('Filtertime (s)')
ylabel('K1')

subplot(2,3,2), hold on
x = K1p(a:z,example_orn);
y = LFP(a:z,example_orn);
plot(x,y,'.')
legend(['r^2 = ' oval(rsquare(x,y))],'Location','southeast')
xlabel('K1 \otimes s(t)')
ylabel('dLFP (mV/s)')
title('No gain correction')

subplot(2,3,5), hold on
x = convolve(1e-3*(1:length(PID)),PID(:,example_orn),Weber_K,ft); x = x(a:z);
y = LFP(a:z,example_orn);
plot(x,y,'.')
legend(['r^2 = ' oval(rsquare(x,y))],'Location','southeast')
xlabel('K1_{Weber} \otimes s(t)')
ylabel('dLFP (mV/s)')
title('No gain correction')

subplot(2,3,3), hold on
S = [K1p(:,example_orn), PID(:,example_orn)];
R = adaptiveGainModel(S,weber_data.p);
plot(R(a:z),LFP(a:z,example_orn),'.')
legend(['r^2 = ' oval(rsquare(R(a:z),LFP(a:z,example_orn)))],'Location','southeast')
xlabel('g_{Weber}(t)(K1 \otimes s(t))')
title('Corrected by Weber-calculated Gain')

subplot(2,3,6), hold on
S = [convolve(1e-3*(1:length(PID)),PID(:,example_orn),Weber_K,ft), PID(:,example_orn)];
R = adaptiveGainModel(S,weber_data.p);
plot(R(a:z),LFP(a:z,example_orn),'.')
legend(['r^2 = ' oval(rsquare(R(a:z),LFP(a:z,example_orn)))],'Location','southeast')
xlabel('g_{Weber}(t)(K1_{Weber} \otimes s(t))')
title('Corrected by Weber-calculated Gain')
ylabel('dLFP (mV/s)')

labelFigure
prettyFig;

if being_published	
	snapnow	
	delete(gcf)
end

%% Naturalistic Stimulus: The LFP to Firing Rate Module
% In this section, we look at the mapping from the LFP to the firing rate. In the following figure, we plot the firing rate of the neuron vs. the filter convolved with the derivative of the LFP (left). We then fit a nonlinearity (a logistic function) to see how much better we can do. 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
x = K2p(:,example_orn); x = x(a:z);
y = fA(:,example_orn); y = y(a:z);
l = plot(x(1:ss:end),y(1:ss:end),'.');
legend(l,['r^2 = ' oval(rsquare(x,y))],'Location','southeast')
xlabel('K2 \otimes LFP')
ylabel('Firing Rate (Hz)')

ft = fittype('logistic(x,A,k,x0)');
ff = fit(x(:),y(:),ft,'StartPoint',[75 .1 -.45]);
plot(sort(x),ff(sort(x)),'r')

subplot(1,2,2), hold on
x = ff(K2p(:,example_orn)); x = x(a:z);
y = fA(:,example_orn); y = y(a:z);
l = plot(x(1:ss:end),y(1:ss:end),'.');
legend(l,['r^2 = ' oval(rsquare(x,y))],'Location','southeast')
xlabel('f(K2 \otimes LFP)')
ylabel('Firing Rate (Hz)')

labelFigure
prettyFig;

if being_published	
	snapnow	
	delete(gcf)
end

%%
% Now, we consider how a dynamic contrast-sensitive term can help with this:
% 

clear d
d.stimulus = [LFP(:,example_orn), K2p(:,example_orn)];
d.response = fA(:,example_orn);
d.response(1:1e4) = NaN;

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
clear p
p. x0 = -0.6189;
p.  A = 74.5937;
p. k0 = 0.1168;
p.  n = 11.9531;
p.tau = 38.8125;
p.  B = 0.0031;

[R,~,K] = contrastLogisticModel(d.stimulus,p);
l = plot(R(a:ss:z),fA(a:ss:z,example_orn),'.');
legend(l,['r^2 = ' oval(rsquare(R(a:z),fA(a:z,example_orn)))],'Location','southeast')
xlabel('f(K2 \otimes dLFP)')
ylabel('Firing Rate (Hz)')

subplot(1,2,2), hold on
title(['\tau_{\sigma} = ' oval(p.tau*p.n), 'ms'])
hx = linspace(0.104,0.118,100);
hy = histcounts(K,hx); hy = hy/sum(hy); hx = hx(1:end-1) + mean(diff(hx));
plot(hx,hy);
xlabel('Steepness parameter')
ylabel('Probability')

labelFigure
prettyFig('FixLogX',true);

if being_published	
	snapnow	
	delete(gcf)
end
 
%%
% How critical is this timescale? In the following figure, we vary the timescale and see how if affects the fit:

all_t = round(logspace(.2,log10(10e3),50));
r2 = NaN*all_t;
for i = 1:length(all_t)
	p.tau = all_t(i);
	R = contrastLogisticModel(d.stimulus,p);
	r2(i) = rsquare(R(a:z),fA(a:z,example_orn));
end


figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(all_t,r2,'k+')
set(gca,'XScale','log')
xlabel('\tau_{contrast gain} (ms)')
ylabel('r^2')

labelFigure
prettyFig('FixLogX',true);

if being_published	
	snapnow	
	delete(gcf)
end

%%
% So it looks like the timescale is very poorly constrained. 


%% Naturalistic Stimulus: Firing Rate
% In this section, we look at the LFP to firing rate transformation, and in particular, if our understanding of contrast adaptation helps improve the prediction of the firing rate. In the following figure, we attempt to go from the stimulus to the firing rate, and at each stage, correct by the mean- or contrast-sensitive changes we previously observed. 

% first show the firing vs. the uncorrected projections
ss = 2;
figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[800 800]); hold on
ft = 1e-3*(1:length(K1)) - .1;


% ---------------------  PID -> LFP -> firing rate ------------------------------
subplot(2,2,1), hold on
x = convolve(1e-3*(1:length(PID)),K1p(:,example_orn),K2(:,example_orn),ft);
x = x(a:z);
y = fA(a:z,example_orn);
l = plot(x(1:ss:end),y(1:ss:end),'.');
legend(l,['r^2 = ' oval(rsquare(x,y))],'Location','southeast')
xlabel('K2 \otimes (K1 \otimes s(t))')
ylabel('Firing Rate (Hz)')

% ---------------------  PID -> LFP*Weber -> firing rate ------------------------------
subplot(2,2,2), hold on
S = [convolve(1e-3*(1:length(PID)),PID(:,example_orn),K1(:,example_orn),ft), PID(:,example_orn)];
weber_corrected_LFP = adaptiveGainModel(S,weber_data.p);
% get a filter from this to the firing rate
K = fitFilter2Data(weber_corrected_LFP(a:z),y,'offset',200);
K = K(101:800);
weber_corrected_fp = convolve(1e-3*(1:length(PID)),weber_corrected_LFP,K,ft);
x = weber_corrected_fp(a:z);
l = plot(x(1:ss:end),y(1:ss:end),'.');
legend(l,['r^2 = ' oval(rsquare(x,y))],'Location','southeast')
xlabel('K2 \otimes (g_{Weber}(t)*K1 \otimes s(t)))')
ylabel('Firing Rate (Hz)')

% --------------------- f( PID -> LFP*Weber -> firing rate) ------------------------------

subplot(2,2,3), hold on
d.stimulus = [LFP(:,example_orn),weber_corrected_fp];

clear p
p. x0 = -0.8824;
p.  A = 73.9977;
p. k0 = 0.0732;
p.  n = 32.4375;
p.tau = 4.9219;
p.  B = 0.0082;


[R,~,K] = contrastLogisticModel(d.stimulus,p);
l = plot(R(a:ss:z),fA(a:ss:z,example_orn),'.');
legend(l,['r^2 = ' oval(rsquare(R(a:z),fA(a:z,example_orn)))],'Location','southeast')
xlabel('f(K2 \otimes (g_{Weber}(t)*K1 \otimes s(t))))')
ylabel('Firing Rate (Hz)')
title(['\tau_{sigma} = ', oval(p.tau*p.n), 'ms'])

% plot the distribution of k
subplot(2,2,4), hold on
x = linspace(min(K)/2,1.1*max(K),100);
hy = histcounts(K,x);
x = x(1:end-1) + mean(diff(x)); hy = hy/sum(hy);
plot(x,hy)
xlabel('Steepness parameter')
ylabel('Probability')

labelFigure
prettyFig('fs',14);

if being_published	
	snapnow	
	delete(gcf)
end


%% Version Info
%
pFooter;


