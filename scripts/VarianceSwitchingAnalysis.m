% VarianceSwitchingAnalysis.m
% analysis of variance switching experiment
% 
% created by Srinivas Gorur-Shandilya at 2:34 , 27 August 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


% this code determines if this function is being called by publish() or not
calling_func = dbstack;
% add homebrew path
path1 = getenv('PATH');
path1 = [path1 ':/usr/local/bin'];
setenv('PATH', path1);

being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
		unix(['tag -a publish-failed ',which(mfilename)]);
		unix(['tag -r published ',which(mfilename)]);
	end
end
tic

%% Variance Switching Analysis
% In this document, we analyse the responses of ORNs to stimuli where we rapidly switch from two ensembles of stimuli, differing only in their means. 


%     ######  ######## #### ##     ## ##     ## ##       ##     ##  ######  
%    ##    ##    ##     ##  ###   ### ##     ## ##       ##     ## ##    ## 
%    ##          ##     ##  #### #### ##     ## ##       ##     ## ##       
%     ######     ##     ##  ## ### ## ##     ## ##       ##     ##  ######  
%          ##    ##     ##  ##     ## ##     ## ##       ##     ##       ## 
%    ##    ##    ##     ##  ##     ## ##     ## ##       ##     ## ##    ## 
%     ######     ##    #### ##     ##  #######  ########  #######   ######  


%% Stimulus
% In the following figure, we show what the stimulus looks like. On the right, we show the distributions of the stimulus in the two cases. For some reason, the distribution when the variance is high is no longer a nice looking Gaussian, even though that was what we had when we tested it. 

path_name = '/local-data/DA-paper/switching/variance/v2/';
[PID, LFP, fA, paradigm, orn] = consolidateData(path_name,1);

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

% make predictions everywhere using the DA model
p.tau_z =  101.2500;
p.tau_y =  25.9219;
p.  n_y =  2;
p.  n_z =  2;
p.    A =  221.5312;
p.    B =  2.7343;
p.    C =  0.0076;
p.   s0 =  -0.3125;

DA_pred = NaN*fA;
for i = 1:width(fA)
	DA_pred(:,i) = DAModelv2(PID(:,i),p);
end

reshaped_DA = DA_pred(global_start:end-1e4-1,1:width(PID));
reshaped_DA = reshape(reshaped_DA,block_length,width(reshaped_DA)*length(reshaped_DA)/block_length);




% also reshape the orn ID
reshaped_orn = repmat(orn,length(global_start:length(PID)-1e4-1)/block_length,1);
reshaped_orn = reshaped_orn(:);

% throw our NaNs globally. so we're throwing out epochs where the data is incomplete
rm_this = isnan(sum(reshaped_LFP));
reshaped_LFP(:,rm_this) = [];
reshaped_PID(:,rm_this) = [];
reshaped_fA(:,rm_this) = [];
reshaped_DA(:,rm_this) = [];
reshaped_orn(rm_this) = [];

% make colour scheme for block analysis
filter_length = 1000;
offset = 200;
all_offsets = [1 3 6 8];
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


%                  ##       ######## ########  
%                  ##       ##       ##     ## 
%                  ##       ##       ##     ## 
%                  ##       ######   ########  
%                  ##       ##       ##        
%                  ##       ##       ##        
%                  ######## ##       ##        


%% LFP Analysis
% We now line up all the LFP signals by the switching time, and look at how the LFP changes in aggregate when we switch from a low variance stimulus to a high variance stimulus. We band pass the LFP to remove spikes and slow fluctuations irrelevant to this analysis. In the following figure, we plot all the traces, together with the mean. We can clearly see that the variance of the LFP follows the variance of the stimulus well, and decreases quickly when we switch to the low variance stimulus. 


% make colour scheme for block analysis
filter_length = 1000;
offset = 200; % for the filter, not the nonlinearity
all_offsets = [0:0.1:10];
window_length = 1; % this is calculated around the offset

time = 1e-3*(1:length(reshaped_LFP));
figure('outerposition',[0 0 600 500],'PaperUnits','points','PaperSize',[1000 500]); hold on

% fake a colourbar
yy0 = mean(mean(reshaped_LFP)) + 30*std(mean2(reshaped_LFP));
c = parula(300);
for i = 1:length(c)
	xx = time([floor(length(time)*(i/length(c))) floor(length(time)*(i/length(c)))]);
	yy = [yy0  yy0*1.1];
	plot(xx,yy,'Color',c(i,:),'LineWidth',4)
end

plot(time,reshaped_LFP(:,1:ss:end),'Color',[.5 .5 .5 .5]);
plot(time,mean2(reshaped_LFP),'Color','k','LineWidth',4);
xlabel('Time since low \rightarrow high switch (s)')
ylabel('\DeltaLFP (mV)')

prettyFig

if being_published
	snapnow
	delete(gcf)
end


%%
% We now extract filters in two second blocks in this triggered time (starting from the time of switch from high to low). In the following figure, we show the filters we extract on the left, with error bars. The middle panel shows how good the filter we extract from each segment is in predicting the LFP. The plot on the right shows the residuals of the data and the linear prediction: changes in slope of this plot correspond to changes in gain at the LFP level. 


% we are going to calculate only one filter/epoch
sr = 1e3; % sampling rate, Hz
if exist('.cache/VSA_K.mat','file') == 2
	load('.cache/VSA_K.mat','K')
else
	filter_length = 1000;
	offset = 200;
	K = NaN(2,filter_length-offset,width(reshaped_LFP));
	for i = 1:width(reshaped_LFP)
		textbar(i,width(reshaped_PID))

		% calculate filter for large variance epoch
		stim = reshaped_PID(:,i);
		resp = reshaped_LFP(:,i);

		resp(1:1e3) = NaN;
		resp(5e3:end)= NaN;

		try
			[this_K,ft] = fitFilter2Data(stim,resp,'reg',1,'offset',offset,'filter_length',filter_length);
			K(1,:,i) = this_K(100:end-101);
		catch 
		end

		% calculate filter for low variance epoch
		stim = reshaped_PID(:,i);
		resp = reshaped_LFP(:,i);

		resp(1:6e3) = NaN;

		try
			[this_K,ft] = fitFilter2Data(stim,resp,'reg',1,'offset',offset,'filter_length',filter_length);
			K(2,:,i) = this_K(100:end-101);
		catch 
		end
	end
	mkdir('.cache')
	save('.cache/VSA_K.mat','K')
end
	
% make linear predictions on the de-trended data using a mean filter averaged over all cases
K_mean = mean(squeeze(mean(K,1)),2);
ft = -99:700;
LFP_pred = NaN*reshaped_LFP;
for i = 1:width(reshaped_LFP)
	LFP_pred(:,i) = convolve(1e-3*(1:length(reshaped_PID)),reshaped_PID(:,i),K_mean,ft);
end

% compute the slopes
if exist('.cache/VSA_LFP_slopes.mat','file') == 2
	load('.cache/VSA_LFP_slopes.mat','n')
else
	n = NaN(width(reshaped_PID),length(all_offsets));
	for i = 1:width(reshaped_PID)
		textbar(i,width(reshaped_LFP))
		for j = 1:length(all_offsets)
			x = LFP_pred(:,i);
			y = reshaped_LFP(:,i);

			a = round(1e3*(all_offsets(j) - window_length/2));
			z = round(1e3*(all_offsets(j) + window_length/2));

			x = circshift(x,length(x) - z );
			y = circshift(y,length(y) - z );

			x = x(end-window_length*1e3:end);
			y = y(end-window_length*1e3:end);

			rm_this = isnan(x) | isnan(y);
			x(rm_this) = [];
			y(rm_this) = [];
			try
				ff = fit(x(:),y(:),'poly1');
				n(i,j) = ff.p1;
			catch
			end
		end
	end
	save('.cache/VSA_LFP_slopes.mat','n')
end


figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1), hold on
filtertime = 1e-3*ft;
clear l
L = {'K_{high}','K_{low}'};

% high and low variance
l(1) = plot(filtertime,mean(squeeze(mean((K(1,:,:)),1)),2),'Color',c(100,:));
l(2) = plot(filtertime,mean(squeeze(mean((K(2,:,:)),1)),2),'Color',c(250,:));
legend(l,L,'Location','southeast')

xlabel('Filter Lag (s)')
ylabel('Filter Amplitude')

subplot(1,3,2), hold on
for j = 1:5:length(all_offsets)
	x = LFP_pred;
	y = reshaped_LFP;

	a = round(1e3*(all_offsets(j) - window_length/2));
	z = round(1e3*(all_offsets(j) + window_length/2));

	x = circshift(x,length(x) - z,1);
	y = circshift(y,length(y) - z,1);

	x = x(end-window_length*1e3:end,:);
	y = y(end-window_length*1e3:end,:);

	ci = max([1 floor(length(c)*(all_offsets(j)*sr)/length(reshaped_LFP))]);
	handles(j) = plotPieceWiseLinear(x(:),y(:),'Color',c(ci,:),'nbins',30);
	delete(handles(j).shade)
	delete(handles(j).line(2:3))
end

xlabel('Projected Stimulus')
ylabel('\DeltaLFP (mV)')

subplot(1,3,3), hold on
% fake a colourbar
yy0 = 3.2;
for i = 1:length(c)
	xx = time([floor(length(time)*(i/length(c))) floor(length(time)*(i/length(c)))]);
	yy = [yy0  yy0*1.01];
	plot(xx,yy,'Color',c(i,:),'LineWidth',4)
end

errorShade(all_offsets,mean(n),sem(n'),'Color','k');

set(gca,'YLim',[2 3.5])
ylabel('LFP Gain (mV/V)')
xlabel('Time since switch (s)')

prettyFig()


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

time = 1e-3*(1:length(reshaped_fA));
figure('outerposition',[0 0 600 500],'PaperUnits','points','PaperSize',[1000 500]); hold on

% fake a colourbar
yy0 = nanmean(nanmean(reshaped_fA)) + 10*nanstd(nanmean(reshaped_fA));
c = parula(300);
for i = 1:length(c)
	xx = time([floor(length(time)*(i/length(c))) floor(length(time)*(i/length(c)))]);
	yy = [yy0  yy0*1.01];
	plot(xx,yy,'Color',c(i,:),'LineWidth',4)
end

plot(1e-3*(1:length(reshaped_fA)),reshaped_fA(:,1:ss:end),'Color',[.5 .5 .5 .5]);
plot(1e-3*(1:length(reshaped_fA)),mean2(reshaped_fA),'Color','k','LineWidth',4);
xlabel('Time since low \rightarrow high switch (s)')
ylabel('Firing Rate (Hz)')

prettyFig

if being_published
	snapnow
	delete(gcf)
end


%%
% We now extract LN models in two second blocks in this triggered time (starting from the time of switch from high to low), just like we did with the LFP. In  the figure below, we plot the filters we extract on the left, the residuals showing the output non-linearities in the middle, and exponents of Hill fits on the right. We can clearly see that the slope of the static non-linearity in the middle is higher when the variance is low. In the panel on the right, we quantify this apparent gain increase. 


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
			[this_K2,ft] = fitFilter2Data(stim,resp,'reg',1,'offset',offset,'filter_length',filter_length);
			K2(1,:,i) = this_K2(100:end-101);
		catch 
		end

		% calculate filter for low variance epoch
		stim = reshaped_PID(:,i);
		resp = reshaped_fA(:,i);

		resp(1:6e3) = NaN;

		try
			[this_K2,ft] = fitFilter2Data(stim,resp,'reg',1,'offset',offset,'filter_length',filter_length);
			K2(2,:,i) = this_K2(100:end-101);
		catch 
		end
	end
	mkdir('.cache')
	save('.cache/VSA_K2.mat','K2')
end


% make linear predictions on the de-trended data using a mean filter averaged over all cases
K2_mean = nanmean(squeeze(nanmean(K2,1)),2);
ft = -99:700;
fA_pred = NaN*reshaped_fA;
for i = 1:width(reshaped_fA)
	fA_pred(:,i) = convolve(1e-3*(1:length(reshaped_PID)),reshaped_PID(:,i),K2_mean,ft);
end

% compute the slopes
if exist('.cache/VSA_fA_slopes.mat','file') == 2
	load('.cache/VSA_fA_slopes.mat','n')
else
	n = NaN(width(reshaped_PID),length(all_offsets));
	for i = 1:width(reshaped_PID)
		textbar(i,width(reshaped_fA))
		for j = 1:length(all_offsets)
			x = fA_pred(:,i);
			y = reshaped_fA(:,i);

			a = round(1e3*(all_offsets(j) - window_length/2));
			z = round(1e3*(all_offsets(j) + window_length/2));

			x = circshift(x,length(x) - z );
			y = circshift(y,length(y) - z );

			x = x(end-window_length*1e3:end);
			y = y(end-window_length*1e3:end);

			rm_this = isnan(x) | isnan(y);
			x(rm_this) = [];
			y(rm_this) = [];

			y = y(x > 0.33*max(x) & x < .66*max(x));
			x = x(x > 0.33*max(x) & x < .66*max(x));

			try
				ff = fit(x(:),y(:),'poly1');
				n(i,j) = ff.p1;
			catch
			end
		end
	end
	save('.cache/VSA_fA_slopes.mat','n')
end


figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1), hold on
filtertime = 1e-3*ft;
clear l
L = {'K_{high}','K_{low}'};

% high and low variance
l(1) = plot(filtertime,nanmean(squeeze(nanmean((K2(1,:,:)),1)),2),'Color',c(100,:));
l(2) = plot(filtertime,nanmean(squeeze(nanmean((K2(2,:,:)),1)),2),'Color',c(250,:));
legend(l,L,'Location','southeast')

xlabel('Filter Lag (s)')
ylabel('Filter Amplitude')

subplot(1,3,2), hold on
for j = 1:5:length(all_offsets)
	x = fA_pred;
	y = reshaped_fA;

	a = round(1e3*(all_offsets(j) - window_length/2));
	z = round(1e3*(all_offsets(j) + window_length/2));

	x = circshift(x,length(x) - z,1);
	y = circshift(y,length(y) - z,1);

	x = x(end-window_length*1e3:end,:);
	y = y(end-window_length*1e3:end,:);

	x = x(:); y = y(:);
	rm_this = isnan(x) | isnan(y);


	ci = max([1 floor(length(c)*(all_offsets(j)*sr)/length(reshaped_LFP))]);
	handles(j) = plotPieceWiseLinear(x(~rm_this),y(~rm_this),'Color',c(ci,:),'nbins',30);
	delete(handles(j).shade)
	delete(handles(j).line(2:3))
end

xlabel('Projected Stimulus')
ylabel('Firing Rate (Hz)')

subplot(1,3,3), hold on
% fake a colourbar
yy0 = 105;
for i = 1:length(c)
	xx = time([floor(length(time)*(i/length(c))) floor(length(time)*(i/length(c)))]);
	yy = [yy0  yy0*1.01];
	plot(xx,yy,'Color',c(i,:),'LineWidth',4)
end

errorShade(all_offsets,nanmean(n),nanstd(n)/sqrt(width(reshaped_LFP)),'Color','k');

set(gca,'YLim',[50 120])
ylabel('Firing Gain (Hz/V)')
xlabel('Time since switch (s)')

prettyFig()


if being_published
	snapnow
	delete(gcf)
end

%% Optimal Coding
% In this section we investigate the idea that the ORN is doing something like optimal coding, i.e., matching its I/O curve to the statistics of the input. 

figure('outerposition',[0 0 600 800],'PaperUnits','points','PaperSize',[1000 800]); hold on
ax(1) = subplot(2,1,1); hold on
ax(2) = subplot(2,1,2); hold on
xlabel(ax(2),'Projected Stimulus')
ylabel(ax(2),'Normalised Response')
ylabel(ax(1),'Stimulus Probability')

% high variance
temp = fA_pred(1e3:5e3,:); temp = temp(:);
[y,x] = hist(temp,30);
y = y/sum(y);
l(1) = plot(ax(1),x,y,'r');
plot(ax(2),x,cumsum(y),'r--');

x = fA_pred(1e3:5e3,:);
y = reshaped_fA(1e3:5e3,:);

x = x(:); y = y(:);
rm_this = isnan(x) | isnan(y);

[~,data] = plotPieceWiseLinear(x(~rm_this),y(~rm_this),'Color',c(ci,:),'nbins',30,'make_plot',false);
data.y = data.y/max(data.y);
plot(ax(2),data.x,data.y,'r')

% low variance
temp = fA_pred(6e3:end,:); temp = temp(:);
[y,x] = hist(temp,30);
y = y/sum(y);
l(2) = plot(ax(1),x,y,'b');
plot(ax(2),x,cumsum(y),'b--');

x = fA_pred(6e3:end,:);
y = reshaped_fA(6e3:end,:);

x = x(:); y = y(:);
rm_this = isnan(x) | isnan(y);

[~,data] = plotPieceWiseLinear(x(~rm_this),y(~rm_this),'Color',c(ci,:),'nbins',30,'make_plot',false);
data.y = data.y/max(data.y);
plot(ax(2),data.x,data.y,'b')

legend(l,{'High Variance','Low Variance'})

% fake some plots for a nice legend
l(1) = plot(ax(2),NaN,NaN,'k--');
l(2) = plot(ax(2),NaN,NaN,'k');
L = {'Prediction','Data'};
legend(l,L,'Location','southeast')

prettyFig()

if being_published
	snapnow
	delete(gcf)
end

%% DA Model Predictions
% In this section, we fit the DA model to the data and see if it can account for the change in gain with variance. 

% compute the slopes
if exist('.cache/VSA_DA_slopes.mat','file') == 2
	load('.cache/VSA_DA_slopes.mat','n')
else
	n = NaN(width(reshaped_PID),length(all_offsets));
	for i = 1:width(reshaped_PID)
		textbar(i,width(reshaped_DA))
		for j = 1:length(all_offsets)
			x = reshaped_DA(:,i);
			y = reshaped_fA(:,i);

			a = round(1e3*(all_offsets(j) - window_length/2));
			z = round(1e3*(all_offsets(j) + window_length/2));

			x = circshift(x,length(x) - z );
			y = circshift(y,length(y) - z );

			x = x(end-window_length*1e3:end);
			y = y(end-window_length*1e3:end);

			rm_this = isnan(x) | isnan(y);
			x(rm_this) = [];
			y(rm_this) = [];

			y = y(x > 0.33*max(x) & x < .66*max(x));
			x = x(x > 0.33*max(x) & x < .66*max(x));

			try
				ff = fit(x(:),y(:),'poly1');
				n(i,j) = ff.p1;
			catch
			end
		end
	end
	save('.cache/VSA_DA_slopes.mat','n')
end

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
for j = 1:5:length(all_offsets)
	x = reshaped_DA;
	y = reshaped_fA;

	a = round(1e3*(all_offsets(j) - window_length/2));
	z = round(1e3*(all_offsets(j) + window_length/2));

	x = circshift(x,length(x) - z,1);
	y = circshift(y,length(y) - z,1);

	x = x(end-window_length*1e3:end,:);
	y = y(end-window_length*1e3:end,:);

	x = x(:); y = y(:);
	rm_this = isnan(x) | isnan(y);


	ci = max([1 floor(length(c)*(all_offsets(j)*sr)/length(reshaped_LFP))]);
	handles(j) = plotPieceWiseLinear(x(~rm_this),y(~rm_this),'Color',c(ci,:),'nbins',30);
	delete(handles(j).shade)
	delete(handles(j).line(2:3))
end

xlabel('DA Prediction')
ylabel('Firing rate (Hz)')

subplot(1,2,2), hold on
% fake a colourbar
yy0 = 2.4;
for i = 1:length(c)
	xx = time([floor(length(time)*(i/length(c))) floor(length(time)*(i/length(c)))]);
	yy = [yy0  yy0*1.01];
	plot(xx,yy,'Color',c(i,:),'LineWidth',4)
end

errorShade(all_offsets,nanmean(n),sem(n'),'Color','k');

ylabel('Gain w.r.t DA Model')
xlabel('Time since switch (s)')

prettyFig()


if being_published
	snapnow
	delete(gcf)
end


%     ####  ######  ####     ######  ########  ######## 
%      ##  ##    ##  ##     ##    ## ##     ## ##       
%      ##  ##        ##     ##       ##     ## ##       
%      ##   ######   ##     ##       ##     ## ######   
%      ##        ##  ##     ##       ##     ## ##       
%      ##  ##    ##  ##     ##    ## ##     ## ##       
%     ####  ######  ####     ######  ########  ##       

%% Spike Train Statistics
% In this section, we look at the statistics of the spike trains, and see if the distribution of the inter-spike intervals is different in the two cases. Here, we see they are different, which means that they do not adapt perfectly to the two cases. In other words, there are differences in the spike train statistics between low and high variance. (This in contrast to what Brenner, Bialek and van Steveninck saw). 

path_name = '/local-data/DA-paper/switching/variance/v2/';
all_spikes = sparse(2000000,15);
allfiles = dir([path_name '*Variance*.mat']);
c = 1;
for i = 1:length(allfiles)
	textbar(i,length(allfiles))
	clear spikes
	load([path_name allfiles(i).name])
	for j = 1:length(spikes)
		if ~isempty(spikes(j).A)
			for k = 1:width(spikes(j).A)
				if sum(spikes(j).A(k,:)) > 10 && max(diff(find(spikes(j).A(k,:))))*1e-4 < 5  && max((find(spikes(j).A(k,:))))*1e-4 > 190
					all_spikes(:,c) = spikes(j).A(k,:);
					c = c + 1;
				end
			end
		end
	end
end

% lose the ends
all_spikes = all_spikes(20e4:end-10e4-1,:);
s = (ControlParadigm(2).Outputs(3,20e4:end-10e4-1));

% compute ISIs everywhere
lo_isi = zeros(250,width(all_spikes));
hi_isi = zeros(250,width(all_spikes));
for i = 1:width(all_spikes)
	isi = 1e-4*diff(find(all_spikes(s == 0,i)));
	isi = isi/mean(isi);
	[y,x] = histcounts(isi,0:1e-2:2.5);
	y = y/sum(y);
	lo_isi(:,i) = y;

	isi = 1e-4*diff(find(all_spikes(s == 1,i)));
	isi = isi/mean(isi);
	[y,x] = histcounts(isi,0:1e-2:2.5);
	y = y/sum(y);
	hi_isi(:,i) = y;
end

figure('outerposition',[0 0 600 500],'PaperUnits','points','PaperSize',[800 500]); hold on
errorShade(x(2:end),mean2(hi_isi),sem(hi_isi),'Color',[1 0 0]);
errorShade(x(2:end),mean2(lo_isi),sem(lo_isi),'Color',[0 0 1]);
xlabel('ISI (units of mean ISI)')
ylabel('Probability')
prettyFig()


if being_published
	snapnow
	delete(gcf)
end


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
if being_published
	unix(['tag -a published ',which(mfilename)]);
	unix(['tag -r publish-failed ',which(mfilename)]);
end
