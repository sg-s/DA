% fig3_variance_adaptation.m
% makes figure: variance adaptation in ORNs
% 
% created by Srinivas Gorur-Shandilya at 7:10 , 03 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

pHeader;

%% Mutual information
% In this document we check if ORN gain control increases mutual information. The first thing we do is to generate some synthetic data, and measure mutual information between it and various time shifted versions of itself. 

x = randn(1,1e4);
offsets = floor(linspace(0,10,10));
mi = NaN*offsets;
for i = 1:length(offsets)
	mi(i) = information(x,circshift(x',-offsets(i))');
end

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(offsets,mi,'k+')
xlabel('Offset')
ylabel('Mutual Information (bits)')
title('MI between Gaussian noise and its time shifted version')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% Now we repeat this, but with correlated Gaussian noise. 

x = filter(ones(100,1),1,randn(1,1e4));
offsets = floor(linspace(0,100,100));
mi = NaN*offsets;
for i = 1:length(offsets)
	mi(i) = information(x,circshift(x',-offsets(i))');
end

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(offsets,mi,'k+')
xlabel('Offset')
ylabel('Mutual Information (bits)')
title('MI between correlated noise and its time shifted version')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% Now, we generate a synthetic stimulus using the same correlated Gaussian noise. Our response is simply a time-shifted version of the stimulus. We repeat this analysis to see if the mutual information can estimate this delay. 

x = filter(ones(100,1),1,randn(1,1e4));
y = circshift(x',50)';
offsets = floor(linspace(0,100,100));
mi = NaN*offsets;
for i = 1:length(offsets)
	mi(i) = information(x,circshift(y',-offsets(i))');
end

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(offsets,mi,'k+')
xlabel('Offset')
ylabel('Mutual Information (bits)')
title('MI between correlated noise and response delayed by 50 samples')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% OK, that works well. Now we operate on the real stimulus. First, we work on the variance switching experiment. 

dm = dataManager;

%      ########  ########  ######  ##     ##    ###    ########  ######## 
%      ##     ## ##       ##    ## ##     ##   ## ##   ##     ## ##       
%      ##     ## ##       ##       ##     ##  ##   ##  ##     ## ##       
%      ########  ######    ######  ######### ##     ## ########  ######   
%      ##   ##   ##             ## ##     ## ######### ##        ##       
%      ##    ##  ##       ##    ## ##     ## ##     ## ##        ##       
%      ##     ## ########  ######  ##     ## ##     ## ##        ######## 
     
%      ########     ###    ########    ###    
%      ##     ##   ## ##      ##      ## ##   
%      ##     ##  ##   ##     ##     ##   ##  
%      ##     ## ##     ##    ##    ##     ## 
%      ##     ## #########    ##    ######### 
%      ##     ## ##     ##    ##    ##     ## 
%      ########  ##     ##    ##    ##     ## 


[PID, LFP, fA, paradigm, orn, fly] = consolidateData(dm.getPath('e30707e8e8ef6c0d832eee31eaa585aa'),1);

global_start = 40e3; % 40 seconds
global_end = length(PID) - 5e3; 
% bandpass to remove spikes and slow fluctuations
% for i = 1:width(LFP)
% 	a = find(~isnan(LFP(:,i)),1,'first');
% 	z = find(~isnan(LFP(:,i)),1,'last');
% 	LFP(a:z,i) = bandPass(LFP(a:z,i),1000,10)*10; % now in mV
% end

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

% ######## #### ########  #### ##    ##  ######   
% ##        ##  ##     ##  ##  ###   ## ##    ##  
% ##        ##  ##     ##  ##  ####  ## ##        
% ######    ##  ########   ##  ## ## ## ##   #### 
% ##        ##  ##   ##    ##  ##  #### ##    ##  
% ##        ##  ##    ##   ##  ##   ### ##    ##  
% ##       #### ##     ## #### ##    ##  ######   



% we are going to calculate only one filter/epoch
sr = 1e3; % sampling rate, Hz
if exist('.cache/VSA_K2.mat','file') == 2
	load('.cache/VSA_K2.mat','K2')
else
	filter_length = 1000;
	offset = 200;
	K2 = NaN(2,filter_length-offset,width(reshaped_fA));
	for i = 1:width(reshaped_fA)
		textbar(i,width(reshaped_PID))

		% calculate filter for large variance epoch
		stim = reshaped_PID(:,i);
		resp = reshaped_fA(:,i);

		resp(1:1e3) = NaN;
		resp(5e3:end)= NaN;

		try
			this_K2 = fitFilter2Data(stim,resp,'reg',1,'offset',offset,'filter_length',filter_length);
			K2(1,:,i) = this_K2(100:end-101);
		catch 
		end

		% calculate filter for low variance epoch
		stim = reshaped_PID(:,i);
		resp = reshaped_fA(:,i);

		resp(1:6e3) = NaN;

		try
			this_K2 = fitFilter2Data(stim,resp,'reg',1,'offset',offset,'filter_length',filter_length);
			K2(2,:,i) = this_K2(100:end-101);
		catch 
		end
	end
	mkdir('.cache')
	save('.cache/VSA_K2.mat','K2')
end

% compute gains per trial on the uncorrected data
lo_gain = NaN(width(reshaped_PID),1);
hi_gain = NaN(width(reshaped_PID),1);
for i = 1:width(reshaped_PID)
	y = reshaped_fA(1e3:4e3,i);
	x = fA_pred(1e3:4e3,i);
	try
		ff = fit(x(:),y(:),'poly1');
		hi_gain(i) = ff.p1;
	catch
	end

	y = reshaped_fA(6e3:9e3,i);
	x = fA_pred(6e3:9e3,i);
	try
		ff = fit(x(:),y(:),'poly1');
		lo_gain(i) = ff.p1;
	catch
	end
end

lo_gain(lo_gain == 0) = NaN;
hi_gain(hi_gain == 0) = NaN;

% make linear predictions on the de-trended data using a mean filter averaged over all cases
K2_mean = nanmean(squeeze(nanmean(K2,1)),2);
ft = -99:700;
fA_pred = NaN*reshaped_fA;
for i = 1:width(reshaped_fA)
	fA_pred(:,i) = convolve(1e-3*(1:length(reshaped_PID)),reshaped_PID(:,i),K2_mean,ft);
end

%%
% First, we compute the mutual information between the stimulus and the time-shifted stimulus in the data. 

offsets = floor(linspace(-200,200,50));
mi_stim = NaN(length(offsets),width(reshaped_PID));
for i = 1:width(reshaped_PID)
	x = reshaped_PID(:,i)';
	for j = 1:length(offsets)
		mi_stim(j,i) = information(x,circshift(x',-offsets(j))');
	end
end

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
errorShade(offsets,mean(mi_stim,2),sem(mi_stim'));
xlabel('Offset (ms)')
ylabel('Mutual Information (bits)')
title('MI between odor stimulus and its time shifted version')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end


%%
% Now, we compute the mutual information between the stimulus and the firing rate, and compare that to the linear prediction. 

offsets = floor(linspace(-100,200,50));
if ~exist('mi_fA','var')
	mi_fA = NaN(length(offsets),width(reshaped_PID));
	mi_fp = NaN(length(offsets),width(reshaped_PID));
	for i = 1:width(reshaped_PID)
		x = reshaped_PID(:,i)';
		y = reshaped_fA(:,i)';
		if sum(y)>0
			for j = 1:length(offsets)
				mi_fA(j,i) = information(x,circshift(y',-offsets(j))');
			end
		end

		y = fA_pred(:,i)';
		if nansum(y)>0
			for j = 1:length(offsets)
				mi_fp(j,i) = information(x,circshift(y',-offsets(j))');
			end
		end
	end
end

clear l 
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
errorShade(offsets,nanmean(mi_fA,2),sem(mi_fA'),'Color','r');
l(1) = plot(NaN,NaN,'r.','MarkerSize',30);
errorShade(offsets,nanmean(mi_fp,2),sem(mi_fp'),'Color','b');
l(2) = plot(NaN,NaN,'b.','MarkerSize',30);
xlabel('Offset (ms)')
ylabel('Mutual Information (bits)')
title('MI between odor stimulus and...')
legend(l,{'Firing rate','Projected Stimulus'})


prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% So the mutual information between the stimulus and the projected stimulus is more than the mutual information between the stimulus and the firing rate, which is weird, but makes sense as the projected stimulus is derived entirely from the stimulus. 


%%
% Now, we compare the mutual information between the stimulus and the projected stimulus with and without gain control, to see if the effect of gain control is to increase the mutual information. Here, we multiply each epoch of each projected stimulus trial by the gain measured in that epoch and trial. 

if ~exist('mi_fp2','var')
	mi_fp2 = NaN(length(offsets),width(reshaped_PID));
	for i = 1:width(reshaped_PID)
		x = reshaped_PID(:,i)';
		y = fA_pred(:,i)';
		if nansum(y)>0 && ~isnan(lo_gain(i)) && ~isnan(hi_gain(i))
			% correct the projection by the gain
			y(1:5e3) = y(1:5e3)*hi_gain(i);
			y(5e3+1:end) = y(5e3+1:end)*lo_gain(i);
			for j = 1:length(offsets)
				mi_fp2(j,i) = information(x,circshift(y',-offsets(j))');
			end
		end
	end
end

clear l 
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
errorShade(offsets,nanmean(mi_fA,2),sem(mi_fA'),'Color','r');
l(1) = plot(NaN,NaN,'r.','MarkerSize',30);
errorShade(offsets,nanmean(mi_fp,2),sem(mi_fp'),'Color','b');
l(2) = plot(NaN,NaN,'b.','MarkerSize',30);
errorShade(offsets,nanmean(mi_fp2,2),sem(mi_fp2'),'Color','k');
l(3) = plot(NaN,NaN,'k.','MarkerSize',30);
xlabel('Offset (ms)')
ylabel('Mutual Information (bits)')
title('MI between odor stimulus and...')
legend(l,{'Firing rate','Projected Stimulus','Proj.Stim. +gain'})


prettyFig();

if being_published
	snapnow
	delete(gcf)
end


%%
% So it actually decreases! To see if this makes sense, let's take an extreme example of a stimulus, response and extreme variance gain control that rescales the response perfectly: 


S = filter(ones(100,1),100,randn(1e4,1));
S(S>0) = 1;
S(S<1) = -1;
S(5e3:end) = S(5e3:end)*.5;
R = S;
Rg = R;
Rg(R>0) = 1;
Rg(R<0) = -1;

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(3,1,1); hold on
plot(S,'k')
title('Stimulus')
subplot(3,1,2); hold on
plot(R,'b')
title('Linear Response (same as stimulus)')
subplot(3,1,3); hold on
plot(Rg,'r')
title('Response + perfect gain control')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% In the above example, the stimulus changes in variance after some time, and the linear response tracks this (in this trivial case, it is the same as the stimulus). FInally, we achieve perfect gain control by getting the response to rescale perfectly (a la Nemenman). 

%%
% Now, we measure mutual information between the stimulus and the linear response, and between the stimulus and the linear response with perfect gain control. 

offsets = sort([0 floor(linspace(-100,100,50))]);
nreps = 20;
mi_R = NaN(length(offsets),nreps);
mi_Rg = NaN(length(offsets),nreps);
for i = 1:nreps
	S = filter(ones(100,1),100,randn(1,1e4));
	S(S>0) = 1;
	S(S<1) = -1;
	S(5e3:end) = S(5e3:end)*.5;
	R = S;
	Rg = R;
	Rg(R>0) = 1;
	Rg(R<0) = -1;

	for j = 1:length(offsets)
		mi_R(j,i) = information(S,circshift(R',-offsets(j))');
		mi_Rg(j,i) = information(S,circshift(Rg',-offsets(j))');
	end
end


clear l 
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
errorShade(offsets,nanmean(mi_R,2),sem(mi_R'),'Color','b');
l(1) = plot(NaN,NaN,'b.','MarkerSize',30);
errorShade(offsets,nanmean(mi_Rg,2),sem(mi_Rg'),'Color','r');
l(2) = plot(NaN,NaN,'r.','MarkerSize',30);
xlabel('Offset (ms)')
ylabel('Mutual Information (bits)')
title('MI between synthetic stimulus and...')
legend(l,{'linear response','lin. resp + gain control'})

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% So the effect of gain control is to decrease mutual information. Which makes sense, as with variance gain control, we throw away information on the scale of the fluctuations, retaining only information on *when* the fluctuations happen. 



%% Version Info
% 

pFooter;


