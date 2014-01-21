% look at all data files, compute the LN model and the linear prediction,
% and look at fast adaptation for each, and compute slope of top 10% to bottom
% 10% for all files. 
%% make the output data structures
HistoryLength.max = 5000; % ms
HistoryLength.min = 20; % ms
% scan different history lengths logarithmically
step=(log(HistoryLength.max)-log(HistoryLength.min))/29;
hl = exp(log(HistoryLength.min):step:log(HistoryLength.max));
% make sure it is divisible by 3
hl = 3*floor(hl/3);

%% input data files
allfiles = dir('/data/random-stim/*rand*.mat');
unique_odours = {'5ol','1but','1o3ol','d2succ','2but','i5ac','2ac','5ol'};

%% these are the outputs of this script
LNSlopeRatios = NaN(length(allfiles),length(hl));
LNFitQualityLow = NaN(length(allfiles),length(hl));
LNFitQualityHigh = NaN(length(allfiles),length(hl));
LinearSlopeRatios = NaN(length(allfiles),length(hl));
LinearFitQualityLow = NaN(length(allfiles),length(hl));
LinearFitQualityHigh = NaN(length(allfiles),length(hl));
LinearFitLow = NaN(length(allfiles),length(hl));
LinearFitLowErr = NaN(length(allfiles),length(hl));
LinearFitHigh = NaN(length(allfiles),length(hl));
LinearFitHighErr = NaN(length(allfiles),length(hl));
LNFitLow = NaN(length(allfiles),length(hl));
LNFitLowErr = NaN(length(allfiles),length(hl));
LNFitHigh = NaN(length(allfiles),length(hl));
LNFitHighErr = NaN(length(allfiles),length(hl));
ShortestPulse = NaN(1,length(allfiles));
LongestPulse = NaN(1,length(allfiles));

for i = 1:length(allfiles)
	disp(allfiles(i).name)
	% for each file, calcualte the slopes, etc.
	output_data=AnalyseLinearPrediction(strcat('/data/random-stim/',allfiles(i).name),0,hl);
	
	% parse and combine
	ShortestPulse(i) = output_data(1).ShortestPulse;
	LongestPulse(i) = output_data(1).LongestPulse;

	LinearFitLow (i,:) =output_data(2).low_slopes;
	LinearFitHigh (i,:) =output_data(2).high_slopes;
	LinearFitLowErr (i,:) =output_data(2).low_slopes_err;
	LinearFitHighErr (i,:) =output_data(2).high_slopes_err;

	LNFitLow (i,:) =output_data(1).low_slopes;
	LNFitHigh (i,:) =output_data(1).high_slopes;
	LNFitLowErr (i,:) =output_data(1).low_slopes_err;
	LNFitHighErr (i,:) =output_data(1).high_slopes_err;


	LNSlopeRatios(i,:) = output_data(1).high_slopes./output_data(1).low_slopes;
	LinearSlopeRatios(i,:) = output_data(2).high_slopes./output_data(2).low_slopes;

	LNFitQualityLow(i,:) = output_data(1).low_gof;
	LNFitQualityHigh(i,:) = output_data(1).high_gof;
	LinearFitQualityLow(i,:) = output_data(2).low_gof;
	LinearFitQualityHigh(i,:) = output_data(2).high_gof;

end

% process slope ratios based on fit quality
fit_thresh=  0.7;
LNSlopeRatios(LNFitQualityHigh<fit_thresh)=NaN;
LNSlopeRatios(LNFitQualityLow<fit_thresh)=NaN;
LinearSlopeRatios(LinearFitQualityHigh<fit_thresh)=NaN;
LinearSlopeRatios(LinearFitQualityLow<fit_thresh)=NaN;

% make a plot of the slope ratios for the linear prediction and the LN model 
figure, hold on, subplot(1,2,1), hold on
temp = cmapping(LNSlopeRatios,'jet',[0.5 1.5]);
imagescnan(temp,[0.5 1.5])
title('LN Model','FontSize',20)
xlabel('Window Length (ms)','FontSize',20)
set(gca,'XTick',1:3:length(hl),'XTickLabel',num2cell(hl(1:3:end)))
set(gca,'YLim',[0.5 length(allfiles)+0.5])
set(gca,'XLim',[0.5 length(hl)+0.5])

subplot(1,2,2), hold on
temp = cmapping(LinearSlopeRatios,'jet',[0.5 1.5]);
imagescnan(temp,[0.5 1.5])
title('Linear Prediction','FontSize',20)
xlabel('Window Length (ms)','FontSize',20)
set(gca,'XTick',1:3:length(hl),'XTickLabel',num2cell(hl(1:3:end)))
set(gca,'YLim',[0.5 length(allfiles)+0.5])
set(gca,'XLim',[0.5 length(hl)+0.5])

suptitle('Ratio of Slopes of top 10% to bottom 10%',20)
colorbar

% print a list of IDs and filenames
for i = 1:length(allfiles)
	s = strkat(mat2str(i),'  ',allfiles(i).name);
	disp(s)
end

% plot some things

subplot(1,2,2),hold on
thesefiles=  [11 17]
for i = thesefiles
	temp = LNFitLow(i,:);
	temp2 = LNFitLowErr(i,:);
	temp(LNFitQualityLow(i,:) < 0.7) = NaN;
	temp2(LNFitQualityLow(i,:) < 0.7) = NaN;
	errorbar(hl,temp,temp2,'b')

	temp = LNFitHigh(i,:);
	temp2 = LNFitHighErr(i,:);
	temp(LNFitQualityHigh(i,:) < 0.7) = NaN;
	temp2(LNFitQualityHigh(i,:) < 0.7) = NaN;
	errorbar(hl,temp,temp2,'r')
end