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
% where f(s) is some input nonlinearity. Here, we model the input nonlinearity using a logistic function, and also add a lag to the stimulus. 

%% Synthetic Data
% First, we consider some synthetic data, which we construct by passing the projected LFP through two different logistic functions with different steepnesses in the two epochs. We then back out filters from this synthetic dataset, and show that we do observe contrast adaptation similar to what we see in the real data. 

XG = K2p;
% correct for contrast
for i = 1:width(K2p)
	XG(1e3:5e3,i) = logistic(XG(1e3:5e3,i),55,14,-1.3);
	XG(5e3:9e3,i) = logistic(XG(5e3:9e3,i),55,25,-1.3);
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


figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[500 500]); hold on
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
p.  k = 17.3485;
p. x0 = -0.0441;
p.  d = 0.1065;
p.  A = 5.1093;
p.lag = 52;

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

prettyFig('fs',16)

if being_published
	snapnow
	delete(gcf)
end

%%
% This doesn't seem to work very well. What if we consider the input to the be LFP convolved with the LFP -> firing rate filter? (We're considering the input to by $ K2 \otimes LFP $) 

clear d
these_trials = 50:10:200;
for i = length(these_trials):-1:1
	d(i).stimulus = Ksyn_p(1:end-1e3,these_trials(i));
	d(i).response = XG(1:end-1e3,these_trials(i));
	d(i).response(1:1e3) = NaN;
end


clear p
p.  k = 21.9188;
p. x0 = -1.2980;
p.  d = 0.1061;
p.  A = 5.2343;
p.lag = -10;

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

prettyFig('fs',16)

if being_published
	snapnow
	delete(gcf)
end

%%
% So it isn't doing anything at all, and all improvements come from the fact that we fit a nonlinearity. 

%% Version Info
%
pFooter;