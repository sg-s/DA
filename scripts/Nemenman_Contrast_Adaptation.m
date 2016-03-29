% Nemenman_Contrast_Adaptation.m
% 
% created by Srinivas Gorur-Shandilya at 4:36 , 21 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

pHeader;

path_name = '/local-data/DA-paper/switching/variance/v2/';
[PID, LFP, fA, ~, orn] = consolidateData(path_name,1);

global_start = 40e3; % 40 seconds
global_end = length(PID) - 5e3; 
% bandpass to remove spikes and slow fluctuations
for i = 1:width(LFP)
	LFP(:,i) = bandPass(LFP(:,i),Inf,30); % remove spikes
	LFP(:,i) = -1e4*filtfilt(ones(10,1),10,[0; diff(LFP(:,i))]);
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


%% Nemenman's Model
% In this document, we see if Nemenman's simply NL-like model can account for our observed contrast adaptation. Nemenman's model is as follows:
% 
% $$ \dot{r}=f(s)-d\cdot r $$
%
% where f(s) is a step function around the mean of the input. Here, and also add a lag to the stimulus. 

%% Synthetic Data
% First, we consider some synthetic data, which we construct by passing the projected LFP through two different logistic functions with different steepnesses in the two epochs. We then back out filters from this synthetic dataset, and show that we do observe contrast adaptation similar to what we see in the real data. 

XG = K2p;
% correct for contrast
for i = 1:width(K2p)
	XG(1:5e3,i) = logistic(XG(1:5e3,i),55,14,-1.3);
	XG(5e3+1:end,i) = logistic(XG(5e3+1:end,i),55,25,-1.3);
end

% back out filters
a = 1e3; z = 9e3;
Ksyn = extractFilters(K1p,XG,'use_cache',true,'a',a,'z',z);

% average filters and project stimulus
Ksyn = nanmean(Ksyn,2);


Ksyn_p = NaN*reshaped_fA;
for i = 1:width(reshaped_fA)
	Ksyn_p(:,i) = convolve(1e-3*(1:length(reshaped_LFP)),K1p(:,i),Ksyn,ft);
end


figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plotPieceWiseLinear(Ksyn_p(1e3:5e3,:),XG(1e3:5e3,:),'nbins',50,'Color','r');
plotPieceWiseLinear(Ksyn_p(6e3:end-200,:),XG(6e3:end-200,:),'nbins',50,'Color','b');
xlabel('K \otimes s')
ylabel('Synthetic Data (Hz)')

prettyFig('fs',16)

if being_published
	snapnow
	delete(gcf)
end

% We now fit a Nemenman Model to this dataset. We consider the input to be LFP, and the output to be the synthetic data that we passed through a contrast-sensitive logistic function. 


clear d
these_trials = 50:10:200;
for i = length(these_trials):-1:1
	d(i).stimulus = K1p(1:end-1e3,these_trials(i));
	d(i).response = XG(1:end-1e3,these_trials(i));
	d(i).response(1:1e3) = NaN;
end

clear p
p.  d = 0.0095;
p.  A = 0.5760;
p.lag = 15;

% make the prediction using the Nemenman Model
Np = K2p;
for i = 1:width(K2p)
	Np(:,i) = NemenmanModel(K1p(:,i),p);
end

figure('outerposition',[0 0 900 800],'PaperUnits','points','PaperSize',[900 800]); hold on
subplot(2,2,1), hold on
plotPieceWiseLinear(Ksyn_p(1e3:5e3,:),XG(1e3:5e3,:),'nbins',50,'Color','r');
plotPieceWiseLinear(Ksyn_p(6e3:end-200,:),XG(6e3:end-200,:),'nbins',50,'Color','b');
xlabel('K \otimes s')
ylabel('Synthetic Data (Hz)')

subplot(2,2,2), hold on
plotPieceWiseLinear(Np(1e3:5e3,:),XG(1e3:5e3,:),'nbins',50,'Color','r');
plotPieceWiseLinear(Np(6e3:end-200,:),XG(6e3:end-200,:),'nbins',50,'Color','b');
xlabel('Nemenman Model Prediction')
ylabel('Synthetic Data (Hz)')

subplot(2,2,3); hold on
plotPieceWiseLinear(Ksyn_p(1e3:5e3,:),Np(1e3:5e3,:),'nbins',50,'Color','r');
plotPieceWiseLinear(Ksyn_p(6e3:end-200,:),Np(6e3:end-200,:),'nbins',50,'Color','b');
ylabel('Nemenman Model Prediction')
xlabel('K \otimes s')

subplot(2,2,4); hold on
r2_X = NaN(width(K2p),1);
r2_XG = r2_X;
for i = 1:width(K2p)
	fp = Ksyn_p([1e3:5e3 6e3:9.7e3],i); r = XG([1e3:5e3 6e3:9.7e3],i);
	try
		r2_X(i) = rsquare(fp,r);
	catch
	end
	fp = Np([1e3:5e3 6e3:9.7e3],i);
	try
		r2_XG(i) = rsquare(fp,r);
	catch
	end
end
plot([0 1],[0 1],'k--')
plot(r2_X,r2_XG,'k+')
xlabel('r^2 Linear prediction')
ylabel('r^2 Nemenman Model')

labelFigure
suptitle('Nemenman Model : stimulus is LFP')

labelFigure
prettyFig('fs',16)

if being_published
	snapnow
	delete(gcf)
end

%%
% The model achieves a sort of contrast adaptation, by working as a "bang-bang" system -- it exponentially relaxes from two extremal values. The following figure shows the model fit vs. the synthetic data, showing this "bang-bang" property. 

i = 166;
figure('outerposition',[0 0 1200 700],'PaperUnits','points','PaperSize',[1200 700]); hold on
time = 1e-3*(1:length(XG));
subplot(2,3,1:2), hold on
plot(time,XG(:,i),'k')
plot(time,Np(:,i),'r')
legend('Synthetic Data',['Nemenman Model r^2=', oval(rsquare(XG(:,i),Np(:,i)))])
xlabel('Time (s)')
ylabel('Response (Hz)')
set(gca,'XLim',[1 9])

subplot(2,3,3), hold on
plot(Np(1e3:9e3,i),XG(1e3:9e3,i),'k.')
xlabel('Nemenman Model')
ylabel('Synthetic Data')

subplot(2,3,4:5), hold on
[ax] = plotyy(time,XG(:,i),time,Ksyn_p(:,i));
legend('Synthetic Data',['Linear Fit r^2=', oval(rsquare(XG(:,i),Ksyn_p(:,i)))])
xlabel('Time (ms)')
ylabel('Response (Hz)')
ax(1).XLim = [1 9];
ax(2).XLim = [1 9];
ax(2).YLim = [-.1 .5];

subplot(2,3,6), hold on
plot(Ksyn_p(1e3:9e3,i),XG(1e3:9e3,i),'k.')
xlabel('Linear projection')
ylabel('Synthetic Data')

prettyFig('fs',16)

if being_published
	snapnow
	delete(gcf)
end

%%
% What if we consider the input to the be LFP convolved with the LFP -> firing rate filter? 

clear d
these_trials = 50:10:200;
for i = length(these_trials):-1:1
	d(i).stimulus = Ksyn_p(1:end-1e3,these_trials(i));
	d(i).response = XG(1:end-1e3,these_trials(i));
	d(i).response(1:1e3) = NaN;
end


clear p
p.  d = 0.0093;
p.  A = 0.5559;
p.lag = -56;

% make the prediction using the Nemenman Model
Np = K2p;
for i = 1:width(K2p)
	Np(:,i) = NemenmanModel(Ksyn_p(:,i),p);
end

figure('outerposition',[0 0 900 800],'PaperUnits','points','PaperSize',[900 800]); hold on
subplot(2,2,1), hold on
plotPieceWiseLinear(Ksyn_p(1e3:5e3,:),XG(1e3:5e3,:),'nbins',50,'Color','r');
plotPieceWiseLinear(Ksyn_p(6e3:end-200,:),XG(6e3:end-200,:),'nbins',50,'Color','b');
xlabel('K \otimes s')
ylabel('Synthetic Data (Hz)')

subplot(2,2,2), hold on
plotPieceWiseLinear(Np(1e3:5e3,:),XG(1e3:5e3,:),'nbins',50,'Color','r');
plotPieceWiseLinear(Np(6e3:end-200,:),XG(6e3:end-200,:),'nbins',50,'Color','b');
xlabel('Nemenman Model Prediction')
ylabel('Synthetic Data (Hz)')

subplot(2,2,3); hold on
plotPieceWiseLinear(Ksyn_p(1e3:5e3,:),Np(1e3:5e3,:),'nbins',50,'Color','r');
plotPieceWiseLinear(Ksyn_p(6e3:end-200,:),Np(6e3:end-200,:),'nbins',50,'Color','b');
ylabel('Nemenman Model Prediction')
xlabel('K \otimes s')

subplot(2,2,4); hold on
r2_X = NaN(width(K2p),1);
r2_XG = r2_X;
for i = 1:width(K2p)
	fp = Ksyn_p([1e3:5e3 6e3:9.7e3],i); r = XG([1e3:5e3 6e3:9.7e3],i);
	try
		r2_X(i) = rsquare(fp,r);
	catch
	end
	fp = Np([1e3:5e3 6e3:9.7e3],i);
	try
		r2_XG(i) = rsquare(fp,r);
	catch
	end
end
plot([0 1],[0 1],'k--')
plot(r2_X,r2_XG,'k+')
xlabel('r^2 Linear prediction')
ylabel('r^2 Nemenman Model')

suptitle('Nemenman Model : stimulus is projected LFP')

labelFigure
prettyFig('fs',16)

if being_published
	snapnow
	delete(gcf)
end



%% Application to Natural Stimulus
% In this section, we see if this model can help us understand the responses to naturalistic stimuli. First, we see if we can use the model to go from the LFP to the firing rate. 


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

clear d
d.stimulus = LFP(:,2);
d.response = fA(:,2);
d.response(1:5e3) = NaN;

clear p
p.  d = 0.0085;
p.  A = 0.6440;
p.lag = 3;

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1:2), hold on
time = 1e-3*(1:length(fA));
plot(time,fA(:,2),'k')
R = NemenmanModel(d.stimulus,p);
plot(time,R,'r')
set(gca,'XLim',[20 30])

subplot(1,3,3), hold on
plot(R(a:10:z),fA(a:10:z,2),'k.')
xlabel('Nemenman Model')
ylabel('ORN Firing Rate (Hz)')
legend(['r^2 = ' oval(rsquare(R(a:10:z),fA(a:10:z,2)))],'Location','southeast')

labelFigure
prettyFig('fs',16)

if being_published
	snapnow
	delete(gcf)
end


%% Version Info
%
pFooter;