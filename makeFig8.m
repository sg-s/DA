% makeFig8.m
% makes figure 8: contrast adaptation in ORNs
% 
% created by Srinivas Gorur-Shandilya at 7:10 , 03 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

% add homebrew path
path1 = getenv('PATH');
if isempty(strfind(path1,[':/usr/local/bin']))
    path1 = [path1 ':/usr/local/bin'];
end
setenv('PATH', path1);

% this code determines if this function is being called 
calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
		unix(['tag -a publish-failed ',which(mfilename)]);
		unix(['tag -r published ',which(mfilename)]);
	end
end
tic


path_name = '/local-data/DA-paper/switching/variance/v2/';
[PID, LFP, fA, paradigm, orn] = consolidateData(path_name,1);

global_start = 40e3; % 40 seconds
global_end = length(PID) - 5e3; 
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


% also reshape the orn ID
reshaped_orn = repmat(orn,length(global_start:length(PID)-1e4-1)/block_length,1);
reshaped_orn = reshaped_orn(:);

% throw our NaNs globally. so we're throwing out epochs where the data is incomplete
rm_this = isnan(sum(reshaped_LFP));
reshaped_LFP(:,rm_this) = [];
reshaped_PID(:,rm_this) = [];
reshaped_fA(:,rm_this) = [];
reshaped_orn(rm_this) = [];

% make colour scheme for block analysis
filter_length = 1000;
offset = 200;

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

all_offsets = [1 3 6 8];
window_length = 1;

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


c = parula(5);
figure('outerposition',[0 0 500 900],'PaperUnits','points','PaperSize',[500 900]); hold on
ax(1) = subplot(5,2,1:2); hold on
for i = 3:10
	ax(i-1) = subplot(5,2,i);
	hold(ax(i-1),'on')
end

% show the stimulus on top
time = 1e-3*(1:length(PID));
time = time(global_start:global_start+40e3);
y = PID(global_start:global_start+40e3,1);
ss = 50;
plot(ax(1),time(1:ss:end),y(1:ss:end),'k')

%                  ##       ######## ########  
%                  ##       ##       ##     ## 
%                  ##       ##       ##     ## 
%                  ##       ######   ########  
%                  ##       ##       ##        
%                  ##       ##       ##        
%                  ######## ##       ##        


% show the reshaped LFP
ss=  10;
time = 1e-3*(1:length(reshaped_LFP));
plot(ax(2),time,reshaped_LFP(:,1:ss:end),'Color',[.5 .5 .5 .5]);
plot(ax(2),time,mean2(reshaped_LFP),'Color','k','LineWidth',4);



% show the LFP filters
filtertime = 1e-3*ft;
clear l
L = {'K_{high}','K_{low}'};

% high and low variance
K_high = mean(squeeze(mean((K(1,:,:)),1)),2);
K_low = mean(squeeze(mean((K(2,:,:)),1)),2);

% normalise correctly 
K_high = K_high*sqrt(1/sum(K_high.^2));
K_low = K_low*sqrt(1/sum(K_low.^2));

l(1) = plot(ax(4),filtertime,K_high,'Color',c(1,:));
l(2) = plot(ax(4),filtertime,K_low,'Color',c(4,:));
legend(l,L,'Location','southeast')


for j = 1:length(all_offsets)
	x = LFP_pred;
	y = reshaped_LFP;

	a = round(1e3*(all_offsets(j) - window_length/2));
	z = round(1e3*(all_offsets(j) + window_length/2));

	x = circshift(x,length(x) - z,1);
	y = circshift(y,length(y) - z,1);

	x = x(end-window_length*1e3:end,:);
	y = y(end-window_length*1e3:end,:);

	x = x(:); y = y(:);
	rm_this = isnan(x) | isnan(y);


	ci = max([1 floor(length(c)*(all_offsets(j)*sr)/length(reshaped_LFP))]);
	[~,data] = plotPieceWiseLinear(x(~rm_this),y(~rm_this),'make_plot',false,'nbins',30);
	plot(ax(6),data.x,data.y,'Color',c(j,:))
end

% show how gain changes with time
all_offsets = [0:0.1:10];
errorShade(ax(8),all_offsets,nanmean(n),nanstd(n)/sqrt(width(reshaped_fA)),'Color','k');



% ######## #### ########  #### ##    ##  ######      ########     ###    ######## ########  ######  
% ##        ##  ##     ##  ##  ###   ## ##    ##     ##     ##   ## ##      ##    ##       ##    ## 
% ##        ##  ##     ##  ##  ####  ## ##           ##     ##  ##   ##     ##    ##       ##       
% ######    ##  ########   ##  ## ## ## ##   ####    ########  ##     ##    ##    ######    ######  
% ##        ##  ##   ##    ##  ##  #### ##    ##     ##   ##   #########    ##    ##             ## 
% ##        ##  ##    ##   ##  ##   ### ##    ##     ##    ##  ##     ##    ##    ##       ##    ## 
% ##       #### ##     ## #### ##    ##  ######      ##     ## ##     ##    ##    ########  ######  



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

all_offsets = [1 3 6 8];
window_length = 2;

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


% show the reshaped fA
ss=  10;
time = 1e-3*(1:length(reshaped_fA));
plot(ax(3),time,reshaped_fA(:,1:ss:end),'Color',[.5 .5 .5 .5]);
plot(ax(3),time,mean2(reshaped_fA),'Color','k','LineWidth',4);


% show the fA filters
filtertime = 1e-3*ft;
clear l
L = {'K_{high}','K_{low}'};

% high and low variance
K_high = nanmean(squeeze(nanmean((K2(1,:,:)),1)),2);
K_low = nanmean(squeeze(nanmean((K2(2,:,:)),1)),2);

% normalise correctly 
K_high = K_high*sqrt(1/sum(K_high.^2));
K_low = K_low*sqrt(1/sum(K_low.^2));

l(1) = plot(ax(5),filtertime,K_high,'Color',c(1,:));
l(2) = plot(ax(5),filtertime,K_low,'Color',c(4,:));
legend(l,L,'Location','northeast')

all_offsets = [1 3 6 8];
for j = 1:length(all_offsets)
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


	ci = max([1 floor(length(c)*(all_offsets(j)*sr)/length(reshaped_fA))]);
	[~,data] = plotPieceWiseLinear(x(~rm_this),y(~rm_this),'make_plot',false,'nbins',30);
	plot(ax(7),data.x,data.y,'Color',c(j,:))
end

% show how gain changes with time
all_offsets = [0:0.1:10];
errorShade(ax(9),all_offsets,nanmean(n),nanstd(n)/sqrt(width(reshaped_fA)),'Color','k');



%       ######   #######   ######  ##     ## ######## ######## ####  ######   ######  
%      ##    ## ##     ## ##    ## ###   ### ##          ##     ##  ##    ## ##    ## 
%      ##       ##     ## ##       #### #### ##          ##     ##  ##       ##       
%      ##       ##     ##  ######  ## ### ## ######      ##     ##  ##        ######  
%      ##       ##     ##       ## ##     ## ##          ##     ##  ##             ## 
%      ##    ## ##     ## ##    ## ##     ## ##          ##     ##  ##    ## ##    ## 
%       ######   #######   ######  ##     ## ########    ##    ####  ######   ######  


% cosmetics
xlabel(ax(1),'Time (s)')
ylabel(ax(1),'Stimulus (V)')

xlabel(ax(2),'Time since switch (s)')
ylabel(ax(2),'\DeltaLFP (mV)')

xlabel(ax(3),'Time since switch (s)')
ylabel(ax(3),'Firing Rate (Hz)')

xlabel(ax(4),'Filter Lag (s)')
ylabel(ax(4),'Filter Amplitude')

xlabel(ax(5),'Filter Lag (s)')
ylabel(ax(5),'Filter Amplitude')

ylabel(ax(6),'\DeltaLFP (mV)')
xlabel(ax(6),'Projected Stimulus')

ylabel(ax(7),'Firing Rate (Hz)')
xlabel(ax(7),'Projected Stimulus')

set(ax(8),'YLim',[2 3.5])
ylabel(ax(8),'LFP Gain (mV/V)')
xlabel(ax(8),'Time since switch (s)')

set(ax(9),'YLim',[50 120])
ylabel(ax(9),'ORN Gain (Hz/V)')
xlabel(ax(9),'Time since switch (s)')


prettyFig('fs=12;','plw=1;','lw=1.5;')

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
