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
[h, data_hi] = plotPieceWiseLinear(x,y,'nbins',50,'Color',[1 0 0]);
delete(h.line(2:3))
delete(h.shade)
x = K2p(6e3:end-200,:);
y = reshaped_fA(6e3:end-200,:);
rm_this = (isnan(sum(y)) | isnan(sum(x)));
x(:,rm_this) = []; y(:,rm_this) = []; 
[h, data_lo] = plotPieceWiseLinear(x,y,'nbins',50,'Color',[0 0 1]);
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

%%
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
	d(i).stimulus = K2p(1:end-1e3,these_trials(i));
	d(i).response = reshaped_fA(1:end-1e3,these_trials(i));
	d(i).response(1:1e3) = NaN;
end

clear p
p. x0 = -0.9570;
p.  A = 50.1992;
p. k0 = 92.3486;
p.  n = 2;
p.tau = 400;
p.  B = 53.5000;



XG = K2p;
for i = 1:width(K2p)
	XG(:,i) = contrastLogisticModel(K2p(:,i),p);
end


figure('outerposition',[0 0 900 800],'PaperUnits','points','PaperSize',[900 800]); hold on
subplot(2,2,1), hold on
plotPieceWiseLinear(K2p(1e3:5e3,:),reshaped_fA(1e3:5e3,:),'nbins',50,'Color','r');
plotPieceWiseLinear(K2p(6e3:end,:),reshaped_fA(6e3:end,:),'nbins',50,'Color','b');
xlabel('K \otimes s')
ylabel('Response (Hz)')

subplot(2,2,2), hold on
plotPieceWiseLinear(XG(1e3:5e3,:),reshaped_fA(1e3:5e3,:),'nbins',50,'Color','r');
plotPieceWiseLinear(XG(6e3:end,:),reshaped_fA(6e3:end,:),'nbins',50,'Color','b');
xlabel('$\hat{R}$','interpreter','latex')
ylabel('Response (Hz)')

subplot(2,2,3); hold on
plotPieceWiseLinear(K2p(1e3:5e3,:),XG(1e3:5e3,:),'nbins',50,'Color','r');
plotPieceWiseLinear(K2p(6e3:end,:),XG(6e3:end,:),'nbins',50,'Color','b');
ylabel('$\hat{R}$','interpreter','latex')
xlabel('K \otimes s')

subplot(2,2,4); hold on
r2_X = NaN(width(K2p),1);
r2_XG = r2_X;
for i = 1:width(K2p)
	fp = K2p([1e3:5e3 6e3:10e3],i); r = reshaped_fA([1e3:5e3 6e3:10e3],i);
	try
		r2_X(i) = rsquare(fp,r);
	catch
	end
	fp = XG([1e3:5e3 6e3:10e3],i);
	try
		r2_XG(i) = rsquare(fp,r);
	catch
	end
end
% convert into remaining variance accounted for
[y,x] = histcounts((r2_XG-r2_X)./(1-r2_X),-1:.02:1);
y = y/sum(y);
y = y*length(y);
x = x(1:end-1) + mean(diff(x));
plot(x,y,'k')
plot([0 0],[0 100],'k--')
set(gca,'XLim',[-1 1],'YLim',[0 10])
xlabel(['Additional Variance' char(10) 'accounted for'])
ylabel('pdf')

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

example_orn = 2;
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

weber_corrected_LFP = R;

if being_published	
	snapnow	
	delete(gcf)
end

%% Naturalistic Stimulus: Firing Rate
% In this section, we look at the LFP to firing rate transformation, and in particular, if our understanding of contrast adaptation helps improve the prediction of the firing rate. In the following figure, we plot the firing rate of the neuron vs. various predictors of the firing rate (the first row). From left to right, they are: 
% 
% # The stimulus -> firing rate filter
% # The stimulus -> LFP -> firing rate mapping (2 filters)
% # The LFP -> firing rate filter 
% # The predicted LFP corrected by the Weber gain control term
%

%%
% The second row shows the same plots, but with best-fit contrast-correcting terms. The titles show the best-fit timescales of contrast gain control. Also shown are the shallowest, mean, and steepest non-linearities fit to the data in red in the first row. 

% first show the firing vs. the uncorrected proejctions
ss = 10;
figure('outerposition',[0 0 1500 800],'PaperUnits','points','PaperSize',[1500 800]); hold on

for i = [1:4 6:8]
	ax(i) = subplot(2,4,i); hold on
end

clear d p

% ---------------------  PID -> firing rate ------------------------------

x = K3p(:,example_orn);
x = x(a:z);
y = fA(a:z,example_orn);
l = plot(ax(1),x(1:ss:end),y(1:ss:end),'.');
legend(l,['r^2 = ' oval(rsquare(x,y))],'Location','southeast')
xlabel(ax(1),'K3 \otimes s(t)')
ylabel(ax(1),'Firing Rate (Hz)')
title(ax(1),'No contrast correction')


% ---------------------  PID -> LFP -> firing rate ------------------------------

x = convolve(1e-3*(1:length(PID)),K1p(:,example_orn),K2(:,example_orn),ft);
x = x(a:z);
l = plot(ax(2),x(1:ss:end),y(1:ss:end),'.');
legend(l,['r^2 = ' oval(rsquare(x,y))],'Location','southeast')
xlabel(ax(2),'K2 \otimes (K1 \otimes s(t))')
ylabel(ax(2),'Firing Rate (Hz)')
title(ax(2),'No contrast correction')

d.response = fA(:,example_orn);
d.response(1:a) = NaN; d.response(z:end) = NaN;
d.stimulus = convolve(1e-3*(1:length(PID)),K1p(:,example_orn),K2(:,example_orn),ft);


p. x0 = -0.5821;
p.  A = 69.9455;
p. k0 = 22.8534;
p.  n = 4.1675;
p.tau = 100;
p.  B = 14.1157;

[R,~,k] = contrastLogisticModel(d.stimulus,p);
l = plot(ax(6),R(a:ss:z),fA(a:ss:z,example_orn),'.');
legend(l,['r^2 = ' oval(rsquare(R(a:z),fA(a:z,example_orn)))],'Location','southeast')
xlabel(ax(6),'K2 \otimes (K1 \otimes s(t))')
ylabel(ax(6),'Firing Rate (Hz)')
title(ax(6),['\tau_{sigma} = ' oval(p.tau*p.n), 'ms'])

% also plot the contrast-modualted curves on the previous plot
all_k = [nanmin(k) nanmean(k) nanmax(k)];
x = linspace(nanmin(d.stimulus),nanmax(d.stimulus),100);
for i = 1:length(all_k)
	plot(ax(2),x,logistic(x,p.A,all_k(i),p.x0),'r');
end

% --------------------- LFP -> firing rate ------------------------------

x = K2p(:,example_orn);
x = x(a:z);
l = plot(ax(3),x(1:ss:end),y(1:ss:end),'.');
legend(l,['r^2 = ' oval(rsquare(x,y))],'Location','southeast')
xlabel(ax(3),'K2 \otimes LFP')
ylabel(ax(3),'Firing Rate (Hz)')
title(ax(3),'No contrast correction')

d.stimulus = K2p(:,example_orn);

p.x0  = -0.3321;
p.   A= 71.5392;
p.  k0= 0.1503;
p.   n= 4.9956;
p. tau= 100;
p.   B= 0.0454;

[R,~,k] = contrastLogisticModel(d.stimulus,p);
l = plot(ax(7),R(a:ss:z),fA(a:ss:z,example_orn),'.');
legend(l,['r^2 = ' oval(rsquare(R(a:z),fA(a:z,example_orn)))],'Location','southeast')
xlabel(ax(7),'K2 \otimes LFP')
ylabel(ax(7),'Firing Rate (Hz)')
title(ax(7),['\tau_{sigma} = ' oval(p.tau*p.n), 'ms'])

% also plot the contrast-modualted curves on the previous plot
all_k = [nanmin(k) nanmean(k) nanmax(k)];
x = linspace(nanmin(d.stimulus),nanmax(d.stimulus),100);
for i = 1:length(all_k)
	plot(ax(3),x,logistic(x,p.A,all_k(i),p.x0),'r');
end


% --------------------- Weber-corrected LFP -> firing rate ------------------------------

x = convolve(1e-3*(1:length(PID)),weber_corrected_LFP,K2(:,example_orn),ft);
x = x(a:z);
l = plot(ax(4),x(1:ss:end),y(1:ss:end),'.');
legend(l,['r^2 = ' oval(rsquare(x,y))],'Location','southeast')
xlabel(ax(4),'K2 \otimes Weber-corrected LFP')
ylabel(ax(4),'Firing Rate (Hz)')
title(ax(4),'No contrast correction')

d.stimulus = convolve(1e-3*(1:length(PID)),weber_corrected_LFP,K2(:,example_orn),ft);

p. x0 = -0.9862;
p.  A = 87.9465;
p. k0 = 0.0935;
p.  n = 1.0034;
p.tau = 50.8196;
p.  B = 0.0337;


[R,~,k] = contrastLogisticModel(d.stimulus,p);
l = plot(ax(8),R(a:ss:z),fA(a:ss:z,example_orn),'.');
legend(l,['r^2 = ' oval(rsquare(R(a:z),fA(a:z,example_orn)))],'Location','southeast')
xlabel(ax(8),'K2 \otimes Weber-corrected LFP')
ylabel(ax(8),'Firing Rate (Hz)')
title(ax(8),['\tau_{sigma} = ', oval(p.tau*p.n), 'ms'])

% also plot the contrast-modualted curves on the previous plot
all_k = [nanmin(k) nanmean(k) nanmax(k)];
x = linspace(nanmin(d.stimulus),nanmax(d.stimulus),100);
for i = 1:length(all_k)
	plot(ax(4),x,logistic(x,p.A,all_k(i),p.x0),'r');
end

prettyFig;

if being_published	
	snapnow	
	delete(gcf)
end


%% Version Info
%
pFooter;


