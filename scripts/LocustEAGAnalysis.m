% LocustEAGAnalysis.m
% analysis of data from Zane Aldworth and Mark Stopfer
% 
% created by Srinivas Gorur-Shandilya at 2:34 , 27 August 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-n

% add homebrew path
path1 = getenv('PATH');
if isempty(strfind(path1,':/usr/local/bin'))
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


%% Locust EAG Analysis
% In this document we analyze some data from Zane Aldworth and Mark Stopfer where EAG measurements are made from locust antenna while stimulating them with binary odour signals. 

%% Example Data
% This is what the data looks like. In the following figure, the PID signals have been low-passed filtered by a 30ms filter, and the EAG signals have been high-pass filtered above a 2 second timescale. There are five trials in this dataset, and the shading indicates the standard error of the mean. 

% load data
load('/local-data/DA-paper/locust/example-data.mat')

% clean up, sub-sample to 1ms
PID = PID1; clear PID1
EAG = EAG1; clear EAG1 

PID = PID(:,1:10:end)';
EAG = EAG(:,1:10:end)';
valve = ODR1(:,1:10:end)';
valve(valve<max(max(valve))/2) = 0;
valve(valve>0) = 1;

% set zero
for i = 1:width(PID)
	PID(:,i) = PID(:,i) - mean(PID(1:300,i));
	EAG(:,i) = EAG(:,i) - mean(EAG(1:300,i));
end

% filter
PID = bandPass(PID,Inf,30);
EAG = bandPass(EAG,2e3,Inf);

t = 1e-3*(1:length(PID));

figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
subplot(2,1,1), hold on
errorShade(t,mean(PID,2),sem(PID'),'Color',[0 0 0]);
ylabel('Stimulus (V)')

subplot(2,1,2), hold on
errorShade(t,mean(EAG,2),sem(EAG'),'Color',[0 0 0]);
ylabel('EAG (mV)')
xlabel('Time (s)')
prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%% 
% We now back out the stimulus to EAG filter using standard methods. The following figure shows the filter extracted from all 5 trials. Shading is the standard error of the mean. The filter is integrating, and is wider than stimulus to ORN filters in the fruit fly. 

offset = 500;
[K, EAG_prediction, gain, gain_err] = extractFilters(PID,EAG,'filter_length',2e3,'filter_offset',offset);
t = 1e-3*(1:length(K)) - 1e-3*offset;

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
errorShade(t,mean(K,2),sem(K'),'Color',[0 0 0]);
xlabel('Filter Lag (s)')
ylabel('Filter')
prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% We now use these filters to project the stimulus and make linear predictions of the response. In the following figure, we compare the linear prediction (red) to the actual response in each trial (black). The top row shows the two time traces overlaid, and the bottom row shows the two traces plotted against each other. We can see that the linear prediction does somewhat OK at predicting the EAG response. 

time = 1e-3*(1:length(PID));
figure('outerposition',[0 0 1400 800],'PaperUnits','points','PaperSize',[1400 800]); hold on
for i = 1:5
	subplot(2,5,i), hold on
	plot(time,EAG(:,i),'k')
	plot(time,EAG_prediction(:,i),'r')
	set(gca,'XLim',[5 15])
	title(['Trial # ' oval(i)])
	if i == 1
		ylabel('EAG/prediction')
		xlabel('Time (s)')
	end

	subplot(2,5,i+5), hold on
	plot(EAG_prediction(1e3:10:end,i),EAG(1e3:10:end,i),'k.')
	if i == 1
		ylabel('EAG')
		xlabel('prediction')
	end
	title(['r^2 = ' oval(rsquare(EAG_prediction(1e3:10:end,i),EAG(1e3:10:end,i)))])
end
prettyFig('fs=20;');

if being_published
	snapnow
	delete(gcf)
end

%% Relationship between gain and stimulus
% We now measure the gain in a 500ms window after every valve onset, and see how that correlates with the mean stimulus in some history length. Estimates of gain are computed by fitting a straight line to a plot of response vs. linear prediction over a 500ms window. Gain estimates are accepted only if linear fits has a r-squared value > 0.8.

%%
% In the following figure, the panel on the left shows the gain evaluated at every valve onset as a function of time, showing that there is some variability in the instantaneous gain. What is the origin of this variability? One hypothesis is that this gain depends on the stimulus in some recent past. To check this, I plot the inst. gain vs. the stimulus in the preceding 250ms. The goodness of fit is measured by a Spearman rank correlation coefficient. 

[ons,offs] = computeOnsOffs(valve(:,1));
inst_gain = NaN(length(ons),width(PID));
mean_stim = NaN(length(ons),width(PID));
gain_r2 = NaN(length(ons),width(PID));

for i = 1:width(PID)
	for j = 5:length(ons)-1
		x = EAG_prediction(ons(j):ons(j)+500,i);
		y = EAG(ons(j):ons(j)+500,i);
		try
			ff = fit(x(:),y(:),'poly1');
			inst_gain(j,i) = ff.p1;
			gain_r2(j,i) = rsquare(x,y);
			mean_stim(j,i) = mean(PID(ons(j)-250:ons(j),i));
		catch
		end
	end
end

inst_gain(gain_r2 < .8) = NaN;

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
errorbar(1e-3*ons,nanmean(inst_gain,2),sem(inst_gain'),'k+')
ylabel('Inst. Gain')
xlabel('Time (s)')

subplot(1,2,2), hold on
[~,data] = plotPieceWiseLinear(mean_stim(~isnan(inst_gain)),inst_gain(~isnan(inst_gain)),'nbins',40,'make_plot',false);
l = errorbar(data.x,data.y,data.ye,'k+');
legend(l,['\rho = ' oval(spear(mean_stim(~isnan(inst_gain)),inst_gain(~isnan(inst_gain))))])
xlabel('Mean Stimulus in preceding 250ms (V)')
ylabel('Inst. Gain')

prettyFig();

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

%%
% This file has the following external dependencies:
showDependencyHash(mfilename);

t = toc;

%% 
% This document was built in: 
disp(strcat(oval(t,3),' seconds.'))

% tag the file as being published 

if being_published
	unix(['tag -a published ',which(mfilename)]);
	unix(['tag -r publish-failed ',which(mfilename)]);
end
