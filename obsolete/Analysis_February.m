%% Dynamical Adaptation in ORNs II
% Do ORNs exhibit fast adaptation to a flickering stimulus? Can a simple dynamical adaptation model predict ORN responses to flickering inputs? Here, I take data from Carlotta's random stimulus experiments and first check how well a simple linear prediction from the stimulus compares to the data. Then, I study how the instantaneous gain of the actual ORN response w.r.t to the predicted ORN response varies with time, and try to find a objective filter from the stimulus to this instantaneous gain to correct for systematic errors in the linear prediction.  

% parameters to tune for figure display
font_size = 20;
font_size2 = 18;
marker_size = 10;
marker_size2 = 20;


%% Summary of results so far
% In the previous analysis (Analysis_January.pdf), I looked at one data file and showed that 1) gain is modulated on a fast timescale (~200ms) and that 2) the PID could not predict the gain well enough to improve the linear prediction. 

%% All Data Included in this analysis
% In this document, all the data I have will be included in this analysis. They are:
disp(ls('/data/random-stim/'))

%% Analysis of fast gain modulation
% For each of these data files, we extract the firing rate of the ORN, and the stimulus (the PID signal), build a linear filter from these two, and then predict the ORN firing rate from the PID. 

% parameters
filter_length = 333;
shift_input = 33;
allfiles = dir('/data/random-stim/*.mat');

% make data structures 
history_lengths = [30 102 150 201 300 402 501 600 801 1002];
hl = history_lengths/3; % history lengths better be divisible by 3!
alldata(1).f = [];
alldata(1).fp = [];
alldata(1).K = [];
alldata(1).PID = [];
alldata(1).time = [];
alldata(1).GainAnalysisData = [];

if exist('alldata.mat') == 2
	load('alldata.mat')
else
	for i = 1:length(allfiles)
		textbar(i,length(allfiles))
		% extract the data
		[PID, time, f] = PrepData3(strcat('/data/random-stim/',allfiles(i).name));
		PID = PID(:);
		time = time(:);
		f = f(:);
		ptrend = fit(time,PID,'Poly1'); 
		PID = PID - (ptrend(time) - mean(ptrend(time)));
		ptrend = fit(time,f,'Poly1'); 
		f = f - (ptrend(time) - mean(ptrend(time)));
		shift_input = 33;
		time = time(1:end-shift_input+1);
		PID = PID(shift_input:end);
		f = f(1:end-shift_input+1);
		f(f<0)=0;

		% compute the filter
		K = FindBestFilter(PID,f);

		% make the prediction
		fp = filter(K,1,PID-mean(PID)) + mean(f);

		% do the gain analysis
		shat = NaN(length(hl),length(PID));
		for j = 1:length(hl)
			shat(j,:) = filtfilt(ones(1,hl(j))/hl(j),1,PID);
			shat(j,1:hl(j)) = NaN;
		end
		[GainAnalysisData] = GainAnalysis(f,fp,PID,shat,history_lengths,hl,filter_length,marker_size,marker_size2,font_size,0);

		% save in workspace
		alldata(i).K = K;
		alldata(i).time = time;
		alldata(i).PID = PID;
		alldata(i).f = f;
		alldata(i).fp = fp;
		alldata(i).GainAnalysisData = GainAnalysisData;
	end
	save('alldata.mat','alldata')
end

time = alldata(1).time;
filtertime = 0:mean(diff(time)):filter_length*mean(diff(time));  % this is the time axis for the filter.
filtertime = filtertime - shift_input*(mean(diff(time))); % correct for shifted input

%%
% The following figure shows all the filters computed from each of the files, with the r-square of the prediction to the data indicated in the title. 
figure('outerposition',[0 0 1200 800],'PaperUnits','points','PaperSize',[1200 800]); hold on
for i = 1:length(alldata)
	autoplot(length(alldata),i); hold on
	plot(filtertime,alldata(i).K/max(alldata(i).K),'LineWidth',2)
	set(gca,'FontSize',font_size2,'box','on','LineWidth',2,'Xlim',[min(filtertime) max(filtertime)])
	title(oval(rsquare(alldata(i).f,alldata(i).fp),2),'FontSize',font_size)
end

%%
% The following figures show all the ORN responses and the predictions for each case overlaid in red. Notice that in most cases, the predicted firing rates are < 0 when the ORN stops firing. 
for i = 1:length(alldata)
	figure('outerposition',[0 0 800 400],'PaperUnits','points','PaperSize',[1000 400]); hold on
	time = alldata(i).time;
	f = alldata(i).f;
	fp = alldata(i).fp;
	plot(time(1000:2000),f(1000:2000),'b','LineWidth',2), hold on
	plot(time(1000:2000),fp(1000:2000),'r','LineWidth',2)
	set(gca,'XLim',[time(1000) time(2000)],'box','on','LineWidth',1.5,'FontSize',font_size2)
	title(strrep(allfiles(i).name,'_','-'),'FontSize',font_size)
	xlabel('Time (s)','FontSize',font_size)
	ylabel('Firing Rate (Hz)','FontSize',font_size)
end

%%
% We now plot the slopes of the prediction to the data for each data set. 
GainRatios = NaN(length(allfiles),length(hl));
GoF = NaN(length(allfiles),length(hl));
for i = 1:length(allfiles)
	GainRatios(i,:) = alldata(i).GainAnalysisData.low_slopes./alldata(i).GainAnalysisData.high_slopes;
	GoF(i,:) = alldata(i).GainAnalysisData.low_gof.*alldata(i).GainAnalysisData.high_gof;
end

hlstring = {'30','102' ,'150','201','300','402','501','600','801','1002'};
figure('outerposition',[0 0 1200 800],'PaperUnits','points','PaperSize',[1200 800]); hold on
subplot(1,2,1), hold on
imagesc(GainRatios), axis tight, colorbar
axis ij
set(gca,'XTick',1:2:10,'XTickLabel',hlstring(1:2:10))
title('Gain of Low Stim/Gain to High Stim','FontSize',font_size2)
xlabel('History Length (ms)','FontSize',font_size)
ylabel('File','FontSize',font_size)
set(gca,'LineWidth',2,'FontSize',font_size)
subplot(1,2,2), hold on
imagesc(GoF), axis tight, colorbar
axis ij
set(gca,'LineWidth',2,'FontSize',font_size)
title('Goodness of fit','FontSize',font_size2)
set(gca,'XTick',1:2:10,'XTickLabel',hlstring(1:2:10))
xlabel('History Length (ms)','FontSize',font_size)
set(gcf,'renderer','zbuffer')

%% Next Steps
% 
% # Prevent predicted firing rates from going < 0 ? 


%% Docs
% This document was generated by MATLAB's _publish_ function. All files needed to generate this document are on a git repository. To clone onto your machine, use:
%
%   git clone https://srinivasgs@bitbucket.com/srinivasgs/da.git
% 
% You will also need a bunch of functions that this depends on, which are also on git. Use
% 
%   git clone https://srinivasgs@bitbucket.com/srinivasgs/core.git
% 
% to get these. 
% 
% Once you have everything, run these commands to generate this document: 
% 
%   options = struct('showCode',false,'format','latex','imageFormat',...
%   'pdf','figureSnapMethod','print','stylesheet','srinivas_latex.xsl');
%
%   publish2('Analysis_January.m',options);
% 

% close all to remove all extraneous figures, since we are publishing to a document anyway
close all
