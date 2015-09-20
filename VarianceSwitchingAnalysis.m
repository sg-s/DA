% VarianceSwitchingAnalysis.m
% analysis of variance switching experiment
% 
% created by Srinivas Gorur-Shandilya at 2:34 , 27 August 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


% this code determines if this function is being called by publish() or not
calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end
tic

%% Variance Switching Analysis
% In this document, we analyse the responses of ORNs to stimuli where we rapidly switch from two ensembles of stimuli, differing only in their means. 

%% Stimulus
% In the following figure, we show what the stimulus looks like. On the right, we show the distributions of the stimulus in the two cases. For some reason, the distribution when the variance is high is no longer a nice looking Gaussian, even though that was what we had when we tested it. 

p = '/local-data/DA-paper/switching/variance/v2/';
[PID, LFP, fA, paradigm, orn] = consolidateData(p,1);

global_start = 40e3; % 40 seconds
global_end = length(PID) - 5e3; 

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,4,1:3), hold on
time = 1e-3*(1:length(PID));
plot(time,PID(:,1),'k');
set(gca,'XLim',[40 100])
xlabel('Time (s)')
ylabel('Stimulus (V)')

subplot(1,4,4), hold on
x = 0:0.01:1.5;
y = NaN(length(x),width(PID));

for i = 1:width(PID)
	ok = repmat([zeros(5e3,1); ones(5e3,1)],length(PID)/1e4,1);
	ok(1:global_start) = 0; % toss the first 40 seconds
	ok(global_end:end) = 0;
	% trim to where we have data
	z = find(isnan(PID(:,i)),1,'first');
	ok(z:end) = 0; ok = logical(ok);
	y(:,i) = hist(PID(ok,i),x);
	y(:,i) = y(:,i)/sum(y(:,i));

end
errorShade(x,nanmean(y,2),sem(y),'Color',[0 0 1],'Shading',0.4);

for i = 1:width(PID)
	ok = repmat([ones(5e3,1); zeros(5e3,1)],length(PID)/1e4,1);
	ok(1:global_start) = 0; % toss the first 40 seconds
	% trim to where we have data
	z = find(isnan(PID(:,i)),1,'first');
	ok(z:end) = 0; ok = logical(ok);
	y(:,i) = hist(PID(ok,i),x);
	y(:,i) = y(:,i)/sum(y(:,i));

end
errorShade(x,nanmean(y,2),sem(y),'Color',[1 0 0],'Shading',0.4);
xlabel('Stimulus (V)')
ylabel('Probability')

prettyFig()

if being_published
	snapnow
	delete(gcf)
end

%%
% To look at this more closely, we trigger the stimulus by epoch, and plot the mean stimulus over all the data. 

% bandpass to remove spikes and slow fluctuations
for i = 1:width(LFP)
	a = find(~isnan(LFP(:,i)),1,'first');
	z = find(~isnan(LFP(:,i)),1,'last');
	LFP(a:z,i) = bandPass(LFP(a:z,i),1000,10)*10; % now in mV
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


% throw our NaNs globally. so we're throwing out epochs where the data is incomplete
rm_this = isnan(sum(reshaped_LFP));
reshaped_LFP(:,rm_this) = [];
reshaped_PID(:,rm_this) = [];
reshaped_fA(:,rm_this) = [];


% make colour scheme for block analysis
filter_length = 1000;
offset = 200;
all_offsets = [1 3 5.8 7.8];
window_length = 2;

figure('outerposition',[0 0 600 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
c = parula(5);
yy = mean(mean(reshaped_PID)) + mean(5*std(reshaped_PID));
for i = 1:length(all_offsets)
	plot([all_offsets(i) all_offsets(i)+window_length],[yy yy],'Color',c(i,:),'LineWidth',10);
end

ss = 10;
plot(1e-3*(1:length(reshaped_PID)),reshaped_PID(:,1:ss:end),'Color',[.5 .5 .5 .5]);
plot(1e-3*(1:length(reshaped_PID)),mean2(reshaped_PID),'Color','k','LineWidth',4);
xlabel('Time since high \rightarrow low switch (s)')
ylabel('Mean Stimulus (V)')

prettyFig

if being_published
	snapnow
	delete(gcf)
end

%%
% It looks as though the mean in the low-variance epoch is a little higher. Is this true? How different are the means of the stimulus in the two epochs? 

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(0.01*randn(width(reshaped_PID),1) + 2, mean(reshaped_PID(5e3:end,:)),'x','Color','b')
errorbar(2,mean(mean(reshaped_PID(5e3:end,:))),std(mean(reshaped_PID(5e3:end,:))),'k')
plot(0.01*randn(width(reshaped_PID),1) + 1, mean(reshaped_PID(1:5e3,:)),'x','Color','r')
errorbar(1,mean(mean(reshaped_PID(1:5e3,:))),std(mean(reshaped_PID(1:5e3,:))),'k')
set(gca,'XLim',[.5 2.5],'XTick',[1 2],'XTickLabel',{'High Variance','Low Variance'},'XMinorTick','off','YLim',[0 .6])

ylabel('Mean Stimulus (V)')
prettyFig
set(gca,'XMinorTick','off')

if being_published
	snapnow
	delete(gcf)
end

%%
% We can see that despite our best efforts, the means are not the same in the two cases. Since the mean is higher in the low variance case, and if this small mean change contributes to gain change, we hypothesise from our earlier results that the gain in the low variance case will be lower than the gain in the high variance case. In the following sections, we see if this is the case or not. 


%% LFP Analysis
% We now line up all the LFP signals by the switching time, and look at how the LFP changes in aggregate when we switch from a low variance stimulus to a high variance stimulus. We band pass the LFP to remove spikes and slow fluctuations irrelevant to this analysis. In the following figure, we plot all the traces, together with the mean. We can clearly see that the variance of the LFP follows the variance of the stimulus well, and decreases quickly when we switch to the low variance stimulus. 


% make colour scheme for block analysis
filter_length = 1000;
offset = 200;
all_offsets = [1 3 5.8 7.8];
window_length = 2;

figure('outerposition',[0 0 600 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
c = parula(5);
yy = mean(mean(reshaped_LFP)) + mean(5*std(reshaped_LFP));
for i = 1:length(all_offsets)
	plot([all_offsets(i) all_offsets(i)+window_length],[yy yy],'Color',c(i,:),'LineWidth',10);
end


plot(1e-3*(1:length(reshaped_LFP)),reshaped_LFP(:,1:ss:end),'Color',[.5 .5 .5 .5]);
plot(1e-3*(1:length(reshaped_LFP)),mean2(reshaped_LFP),'Color','k','LineWidth',4);
xlabel('Time since high \rightarrow low switch (s)')
ylabel('\DeltaLFP (mV)')

prettyFig

if being_published
	snapnow
	delete(gcf)
end


%%
% We now extract filters in two second blocks in this triggered time (starting from the time of switch from high to low). In the following figure, we show the filters we extract on the left, with error bars. The middle panel shows how good the filter we extract from each segment is in predicting the LFP. The plot on the right shows the residuals of the data and the linear prediction: changes in slope of this plot correspond to changes in gain at the LFP level. 

% let's try to pull out filters from every epoch
sr = 1e3; % sampling rate, Hz
K = cache(dataHash([reshaped_LFP,reshaped_PID]));
ft = -99:700;
if isempty(K)
	K = NaN(length(all_offsets),filter_length-offset,width(reshaped_LFP));
	for i = 1:width(reshaped_LFP)
		textbar(i,width(reshaped_PID))
		for j = 1:length(all_offsets)

			x = NaN(length(reshaped_PID),1);
			y = x;

			a = 1 + all_offsets(j)*sr - filter_length + offset;
			z = all_offsets(j)*sr + offset + window_length*sr;
			this_stim = reshaped_PID(a:z,i);
			% ff = fit((1:length(this_stim))',this_stim,'poly1');
			% this_stim = this_stim - ff(1:length(this_stim));
			x(a:z) = this_stim;

			a = all_offsets(j)*sr;
			z = a + window_length*sr;
			this_resp = reshaped_LFP(a:z,i);
			% ff = fit((1:length(this_resp))',this_resp,'poly1');
			% this_resp = this_resp - ff(1:length(this_resp));
			y(a:z) = this_resp;

			otp = false(length(y),1); % only these points
			otp(a:z) = true;

			try
				[this_K,ft] = fitFilter2Data(x,y,'reg',1,'offset',offset,'OnlyThesePoints',otp,'filter_length',filter_length);
				K(j,:,i) = this_K(100:end-102);
				ft = ft(100:end-102);
			catch err
				disp(i)
				disp(err)
			end
		end
	end
	cache(dataHash([reshaped_LFP,reshaped_PID]),K);
end


% make linear predictions on the detrended data
LFP_pred = NaN*reshaped_LFP;
for i = 1:width(reshaped_LFP)
	for j = 1:length(all_offsets)
		a = all_offsets(j)*1e3;
		z = a + window_length*1e3;
		x = reshaped_PID(:,i);
		this_pred = convolve(1e-3*(1:length(x)),x,squeeze(K(j,:,i)),ft);
		LFP_pred(a:z,i) = this_pred(a:z)  - mean(this_pred(a:z));
	end
end


figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1), hold on
filtertime = 1e-3*ft;
for i = 1:length(all_offsets)
	errorShade(filtertime,mean2(squeeze(K(i,:,:))),sem(squeeze(K(i,:,:))),'Color',c(i,:));
end
xlabel('Filter Lag (s)')
ylabel('Filter Amplitude')

subplot(1,3,2), hold on
for j = 1:length(all_offsets)
	a = all_offsets(j)*sr;
	z = window_length*1e3 + a;

	x = LFP_pred(a:z,:);
	y = reshaped_LFP(a:z,:);
	x = x(:); y = y(:);
	handles(j) = plotPieceWiseLinear(x,y,'Color',c(j,:),'nbins',30);
	delete(handles(j).shade)
	delete(handles(j).line(2:3))
end

xlabel('Linear Prediction')
ylabel('\DeltaLFP (mV)')

subplot(1,3,3), hold on
n = NaN(width(reshaped_PID),length(all_offsets));
for i = 1:width(reshaped_PID)
	for j = 1:length(all_offsets)
		a = all_offsets(j)*sr;
		z = window_length*1e3 + a;
		x = LFP_pred(a:z,i);
		y = reshaped_LFP(a:z,i);
		try
			ff = fit(x(:),y(:),'poly1');
			n(i,j) = ff.p1;
		catch
		end
	end
end

clear L
for i = 1:length(all_offsets)
	errorbar(i,mean2(n(:,i)),sem(n(:,i)),'Color',c(i,:),'LineWidth',4)
	plot(.05*randn(length(n(:,i)),1) + i*ones(length(n(:,i)),1),n(:,i),'x','Color',c(i,:))
	L{i} = [oval(all_offsets(i)) '-' oval(all_offsets(i)+window_length) 's'];
end
ylabel('LFP Gain (mV/V)')
prettyFig()
set(gca,'XTick',[1:length(all_offsets)],'XTickLabel',L,'XLim',[0.7 4.1],'XMinorTick','off')
set(gca,'XTickLabelRotation',45)

if being_published
	snapnow
	delete(gcf)
end

% ######## #### ########  #### ##    ##  ######      ########     ###    ######## ########  ######  
% ##        ##  ##     ##  ##  ###   ## ##    ##     ##     ##   ## ##      ##    ##       ##    ## 
% ##        ##  ##     ##  ##  ####  ## ##           ##     ##  ##   ##     ##    ##       ##       
% ######    ##  ########   ##  ## ## ## ##   ####    ########  ##     ##    ##    ######    ######  
% ##        ##  ##   ##    ##  ##  #### ##    ##     ##   ##   #########    ##    ##             ## 
% ##        ##  ##    ##   ##  ##   ### ##    ##     ##    ##  ##     ##    ##    ##       ##    ## 
% ##       #### ##     ## #### ##    ##  ######      ##     ## ##     ##    ##    ########  ######  


%% Firing Rate Analysis
% We now perform a similar analysis, but for the firing rates. First, we trigger all the firing rates by the low-high variance switch, and look at the data. In the following figure, all the traces are plotted in grey, with the mean superimposed in black. 

figure('outerposition',[0 0 600 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
c = parula(5);
yy = nanmean(nanmean(reshaped_fA)) + nanmean(4*nanstd(reshaped_fA));
for i = 1:length(all_offsets)
	plot([all_offsets(i) all_offsets(i)+window_length],[yy yy],'Color',c(i,:),'LineWidth',10);
end

plot(1e-3*(1:length(reshaped_fA)),reshaped_fA(:,1:ss:end),'Color',[.5 .5 .5 .5]);
plot(1e-3*(1:length(reshaped_fA)),mean2(reshaped_fA),'Color','k','LineWidth',4);
xlabel('Time since high \rightarrow low switch (s)')
ylabel('Firing Rate (Hz)')

prettyFig

if being_published
	snapnow
	delete(gcf)
end


%%
% We now extract LN models in two second blocks in this triggered time (starting from the time of switch from high to low), just like we did with the LFP. In  the figure below, we plot the filters we extract on the left, the residuals showing the output non-linearities in the middle, and exponents of Hill fits on the right. We can clearly see that the slope of the static non-linearity in the middle is higher when the variance is low. In the panel on the right, we quantify this apparent gain increase. 

% let's try to pull out filters from every epoch
sr = 1e3; % sampling rate, Hz
K2 = cache(dataHash([reshaped_fA,reshaped_PID]));
ft = -99:700;
if isempty(K2)
	K2 = NaN(length(all_offsets),filter_length-offset,width(reshaped_fA));
	for i = 1:width(reshaped_fA)
		textbar(i,width(reshaped_PID))
		for j = 1:length(all_offsets)

			x = NaN(length(reshaped_PID),1);
			y = x;

			a = 1 + all_offsets(j)*sr - filter_length + offset;
			z = all_offsets(j)*sr + offset + window_length*sr;
			x(a:z) = reshaped_PID(a:z,i);

			a = all_offsets(j)*sr;
			z = a + window_length*sr;
			y(a:z) = reshaped_fA(a:z,i);

			otp = false(length(y),1); % only these points
			otp(a:z) = true;

			try
				[this_K2,ft] = fitFilter2Data(x,y,'reg',1,'offset',offset,'OnlyThesePoints',otp,'filter_length',filter_length);
				K2(j,:,i) = this_K2(100:end-102);
				ft = ft(100:end-102);
			catch err
				disp(i)
				disp(err)
			end
		end
	end
	cache(dataHash([reshaped_fA,reshaped_PID]),[]);
	cache(dataHash([reshaped_fA,reshaped_PID]),K2);
end


% make linear predictions on the detrended data
fp = NaN*reshaped_fA;
for i = 1:width(reshaped_fA)
	for j = 1:length(all_offsets)
		a = all_offsets(j)*1e3;
		z = a + window_length*1e3;
		x = reshaped_PID(:,i);
		this_pred = convolve(1e-3*(1:length(x)),x,squeeze(K2(j,:,i)),ft);
		fp(a:z,i) = this_pred(a:z)  - mean(this_pred(a:z));
	end
end


figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1), hold on
filtertime = 1e-3*ft;
for i = 1:length(all_offsets)
	errorShade(filtertime,mean2(squeeze(K2(i,:,:))),sem(squeeze(K2(i,:,:))),'Color',c(i,:));
end
xlabel('Filter Lag (s)')
ylabel('Filter Amplitude')

subplot(1,3,2), hold on
for j = 1:length(all_offsets)
	a = all_offsets(j)*sr;
	z = window_length*1e3 + a;

	x = fp(a:z,:);
	y = reshaped_fA(a:z,:);
	x = x(:); y = y(:);
	handles(j) = plotPieceWiseLinear(x,y,'Color',c(j,:),'nbins',30);
	delete(handles(j).shade)
	delete(handles(j).line(2:3))
end

xlabel('Linear Prediction')
ylabel('Firing rate (Hz)')

subplot(1,3,3), hold on
% seed_p.       A= 90.7124;
% seed_p.       k= 0.9428;
% seed_p.y_offset= -1.3130;
% seed_p.x_offset= 0.6865;
% seed_p.       n= 3.1016;

% nbins = 10;
% n = zeros(length(all_offsets),nbins);
% block_size = floor(width(x)/nbins);

% for j = 1:length(all_offsets)
% 	a = all_offsets(j)*sr;
% 	z = window_length*1e3 + a;
% 	x = fp(a:z,:);
% 	y = reshaped_fA(a:z,:);
% 	for i = 1:nbins
% 		disp([j i])
% 		aa = block_size*(i-1) + 1;
% 		zz = aa + block_size;
% 		d.stimulus = x(:,aa:zz);
% 		d.response = y(:,aa:zz);
% 		p = fitModel2Data(@hill5,d,'nsteps',50,'UseParallel',false,'Display','none','make_plot',false,'p0',seed_p);
% 		n(j,i) = p.n;
% 	end
% end

n= cache(dataHash([reshaped_fA fp]));
clear L
for i = 1:length(all_offsets)
	errorbar(i,mean(n(i,:)),std(n(i,:))/2,'Color',c(i,:),'LineWidth',4)
	plot(.05*randn(length(n(i,:)),1) + i*ones(length(n(i,:)),1),n(i,:),'x','Color',c(i,:))
	L{i} = [oval(all_offsets(i)) '-' oval(all_offsets(i)+window_length) 's'];
end
ylabel('Hill exponent')
prettyFig()
set(gca,'XTick',[1:length(all_offsets)],'XTickLabel',L,'XLim',[0.7 4.1],'XMinorTick','off')
set(gca,'XTickLabelRotation',45)

if being_published
	snapnow
	delete(gcf)
end

%%
% Here we show that by switching between a low- and a high-variance stimulus, the neuron appears to increase its firing gain during the low-variance stimulus, even though we don't see a similar increase in gain in the LFP. 


%% Version Info
% The file that generated this document is called:
disp(mfilename)

%%
% and its md5 hash is:
Opt.Input = 'file';
disp(dataHash(strcat(mfilename,'.m'),Opt))

%%
% This file should be in this commit:
[status,m]=unix('git rev-parse HEAD');
if ~status
	disp(m)
end

t = toc;

%% 
% This document was built in: 
disp(strcat(oval(t,3),' seconds.'))

% tag the file as being published 
% add homebrew path
path1 = getenv('PATH');
path1 = [path1 ':/usr/local/bin'];
setenv('PATH', path1);

if being_published
	unix(strjoin({'tag -a published',which(mfilename)}));
end
