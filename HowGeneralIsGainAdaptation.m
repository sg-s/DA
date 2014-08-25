% HowGeneralIsGainAdaptation.m
% summarises gain adaptation in all of carlotta's data. makes plots organised by receptor, odor, correlation length, etc. 
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


% the following section pre-computes the processor-heavy parts and caches them for later use by publish()
load('/local-data/DA-paper/data.mat')
do_these = 2:21;

if ~exist('HowGeneralIsGainAdaptation.mat','file')
	
	n = length(do_these);

	% initialise history lengths to run the analysis on.
	history_lengths=[0:0.03:0.3 0.36:0.06:1 1.2:1.2:5];

	% initialise a matrix for all the linear filters
	Filters = NaN(200,n);

	% initialise a matrix for the parameters of all the non-linear functions
	HillFit = NaN(3,n);

	% initialise a matrix that stores all the slopes (and p-values) we calculate from the Gain Analysis
	low_slopes  = NaN(length(history_lengths),n);
	high_slopes = NaN(length(history_lengths),n);
	low_gof  = NaN(length(history_lengths),n);
	high_gof = NaN(length(history_lengths),n);
	p_values = NaN(length(history_lengths),n);

	% initialise a matrix that stores the r-square of the LN fit
	LNFitQuality = NaN(1,n);

	for i = 1:n

		td = do_these(i);

		% fit a linear filter to data
		[K,~,filtertime] = FindBestFilter(data(td).PID,data(td).ORN,[],'filter_length=199;');
		Filters(:,i) = K;

		LinearFit = mean(data(td).ORN)+convolve(data(td).time,data(td).PID,K,data(td).filtertime);
		LinearFit(LinearFit<0)=0;


		% fit a Hill function post-hoc
		xdata = LinearFit;
		ydata = data(td).ORN;

		% crop it to lose NaNs
		ydata(isnan(xdata)) = [];
		xdata(isnan(xdata)) = [];

		xdata = xdata(:);
		ydata = ydata(:);

		fo=optimset('MaxFunEvals',2000,'Display','none');
		x = lsqcurvefit(@hill,[50 2 2],xdata,ydata,[],[],fo);
		LNFit = hill(x,LinearFit);
		HillFit(:,i) = x;

		LNFitQuality(i) = rsquare(LNFit,data(td).ORN);
		if LNFitQuality(i) < 0.8
			disp('Poor fit')
			beep
			keyboard
		end

		% perform gain analysis
		s = 300; % when we start for the gain analysis
		z = length(data(td).ORN) - 33; % where we end 
		clear x
		x.response = data(td).ORN(s:z);
		x.prediction = LNFit(s:z);
		x.stimulus = data(td).PID(s:z);
		x.time = data(td).time(s:z);
		x.filter_length = 200;

		[p_values(:,i),low_slopes(:,i),high_slopes(:,i),low_gof(:,i),high_gof(:,i)] = GainAnalysis3(x,history_lengths);
		
	end

	% cache locally for use later
	filtertime = filtertime*mean(diff(data(td).time));
	save('HowGeneralIsGainAdaptation.mat','Filters','HillFit','LNFitQuality','high_slopes','low_slopes','p_values','filtertime','history_lengths','low_gof','high_gof')
else
	load('HowGeneralIsGainAdaptation.mat')
end

% corrects for points where the p-value is reported as 0 by the bootstrap. it can't be, it's just too small to measure. 
p_values(p_values==0) = min(nonzeros(p_values(:)));

%% How General is Gain Adaptation?
% This document summaries a gain analysis of all of the binary flickering stimuli we have from Carlotta's experiments. 


                                                                           
 %       ###### # #      ##### ###### #####      ####  #    #   ##   #####  ###### 
 %       #      # #        #   #      #    #    #      #    #  #  #  #    # #      
 %       #####  # #        #   #####  #    #     ####  ###### #    # #    # #####  
 %       #      # #        #   #      #####          # #    # ###### #####  #      
 %       #      # #        #   #      #   #     #    # #    # #    # #      #      
 %       #      # ######   #   ###### #    #     ####  #    # #    # #      ###### 
                                                                           

%% Filter Shape
% The following figure shows the best fit linear filters that we backed out of each data set. All filters are normalised by the peak to compare their shapes. 

NormFilters = Filters;

for i = 1:width(Filters)
	NormFilters(:,i) = Filters(:,i)/max(Filters(:,i));
end

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(filtertime,NormFilters)
xlabel('Lag (s)')
ylabel('Filter Amplitude (norm)')
set(gca,'XLim',[min(filtertime)-0.01 max(filtertime)+0.01])
set(gca,'YLim',[-0.9 1.2])
title('Variance in filter shape')

PrettyFig;

snapnow;
delete(gcf);

%       ##     ## #### ##       ##          ######## ##     ## ##    ##  ######  
%       ##     ##  ##  ##       ##          ##       ##     ## ###   ## ##    ## 
%       ##     ##  ##  ##       ##          ##       ##     ## ####  ## ##       
%       #########  ##  ##       ##          ######   ##     ## ## ## ## ##       
%       ##     ##  ##  ##       ##          ##       ##     ## ##  #### ##       
%       ##     ##  ##  ##       ##          ##       ##     ## ##   ### ##    ## 
%       ##     ## #### ######## ########    ##        #######  ##    ##  ######  

%% Variance of Output Non-linearities 
% In each data set, after backing out the filters, we then fit a 3-parameter hill function to the residuals. The following function shows the shape of all the Hill functions. Each curve is the best fit for one data set. 

x = 0:250;
y = NaN(251,length(HillFit));

for i = 1:length(HillFit)
	y(:,i) = hill(HillFit(:,i),x);
end

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(x,y)
xlabel('Linear Output (Hz)')
ylabel('Nonlinearity Output (Hz)')
set(gca,'XLim',[0 255])
title('Variance in Nonlinearity shape')

PrettyFig;

snapnow;
delete(gcf);

%        ###    ##       ##          ########     ###    ########    ###    
%       ## ##   ##       ##          ##     ##   ## ##      ##      ## ##   
%      ##   ##  ##       ##          ##     ##  ##   ##     ##     ##   ##  
%     ##     ## ##       ##          ##     ## ##     ##    ##    ##     ## 
%     ######### ##       ##          ##     ## #########    ##    ######### 
%     ##     ## ##       ##          ##     ## ##     ##    ##    ##     ## 
%     ##     ## ######## ########    ########  ##     ##    ##    ##     ## 

%% Summary of Gain Analysis in all data
% The following figure shows the results of performing a gain analysis using the LN model's prediction for all the data. In the following figure, each pair of red/green curve comes from one data set. Dots indicate points where the slopes are significantly different. 

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(history_lengths,low_slopes,'Color',[0.5 1 0.5])
plot(history_lengths,high_slopes,'Color',[1 0.5 0.5])
set(gca,'XScale','log')

% now plot the dots where significant
for i = 1:length(HillFit)
	sig = p_values(:,i);
	sig = (sig<0.05);

	scatter(history_lengths(sig),low_slopes(sig,i),500,'g.')
	scatter(history_lengths(sig),high_slopes(sig,i),500,'r.')
end

xlabel('History Length (s)')
ylabel('Relative Gain')

PrettyFig;
snapnow;
delete(gcf);

%%
% However, in some cases, the fits to the clouds of points in the gain analysis may not be very good. The following plot shows the distribution of r-square values of the fits that are used to determine the gain in each of these cases, for the entire dataset:

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
[y,x]=hist(low_gof(:),30);
plot(x,y,'g')
[y,x]=hist(high_gof(:),30);
plot(x,y,'r')
legend({'Low Slopes','High Slopes'})
xlabel('r-square of fit')
ylabel('Count')
title('Some fits are very poor')

PrettyFig;

snapnow;
delete(gcf);

%%
% If we retain only the points where the r-square of the fit is >0.8, we end up retaining the following percent of the low slopes:

disp(100*length(low_gof(low_gof>0.8))/(length(low_gof(:))))

%%
% and the following percent of the high-slopes data:

disp(100*length(high_gof(high_gof>0.8))/(length(high_gof(:))))

%%
% In the following plot, we only retain this data:

low_slopes(low_gof<0.8) = NaN;
high_slopes(high_gof<0.8) = NaN;

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(history_lengths,low_slopes,'.-','Color',[0.5 1 0.5])
plot(history_lengths,high_slopes,'.-','Color',[1 0.5 0.5])
set(gca,'XScale','log')

xlabel('History Length (s)')
ylabel('Relative Gain')

% now plot the dots where significant
for i = 1:length(HillFit)
	sig = p_values(:,i);
	sig = (sig<0.05);

	scatter(history_lengths(sig),low_slopes(sig,i),500,'g.')
	scatter(history_lengths(sig),high_slopes(sig,i),500,'r.')
end

PrettyFig;

snapnow;
delete(gcf);

%%
% In general, are the green slopes (gain following low stimuli) signifcanlty higher than the red slopes (gain follwing high stimuli)? To determine this, we plot, for each pair of points corresponding to a single history length and a single dataset, the difference between the low slopes and the high slopes _vs._ the p-value of the difference between that pair (Bonferroni corrected). 

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
d = low_slopes(:) - high_slopes(:);
plot(p_values(:),d,'k.') 
set(gca,'XScale','log')

% plot a reference line for equal gain
plot([.01 1],[0 0 ],'k--')

% plot a reference line for p = 0.05
plot([.05 .05],[min(d)*2 max(d)*2 ],'k--')

set(gca,'XLim',[min(p_values(:))/2 1.1],'YLim',[-0.5 0.6])


ylabel('Low Slopes - High Slopes')
xlabel('p-value (Bonferroni corr.)')

PrettyFig;

snapnow;
delete(gcf);

%% 
% The reference horizontal line indicates equal slope (equal gain for low and high stimuli) and the reference vertical line indicates a p-value of _p=0.05_. Of the statistically significant differences in gain (left half-plane), there are far more points in the top-left quadrant (gain enhancement to low stimuli) than there are in the bottom-left quadrant (gain suppression to low stimuli). 

p_values(high_gof<0.8) = Inf;
p_values(low_gof<0.8) = Inf;
lowest_p = min(p_values);

%%
% Of all the data analysed, only this dataset did not show significant gain control:

disp(data(1+find(lowest_p>0.05)).original_name)


%     ########  ######## ########          ######  ##     ## ########  ######  ##    ## 
%     ##     ## ##       ##     ##        ##    ## ##     ## ##       ##    ## ##   ##  
%     ##     ## ##       ##     ##        ##       ##     ## ##       ##       ##  ##   
%     ########  ######   ########         ##       ######### ######   ##       #####    
%     ##   ##   ##       ##               ##       ##     ## ##       ##       ##  ##   
%     ##    ##  ##       ##        ###    ##    ## ##     ## ##       ##    ## ##   ##  
%     ##     ## ######## ##        ###     ######  ##     ## ########  ######  ##    ## 

%% Reliability and Reproducibility 
% How reprdocicable is the data in this dataset and how reliable is this analysis?  To check this, we compare the following six pieces of data, obtained on different dates with different neurons: 

plothese = [8 10 14 16 17];
for i = 1:length(plothese)
	disp(data(plothese(i)).original_name)
end

%%
% but which all use the same odor, presented in the same pattern, and the neuron measured from is of the same type (ab3A). In the following figure, we plot all the stimuli recorded on all these differnet days, plotted on on top of another: 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
allstim = zeros(length(data(plothese(1)).PID),length(plothese));
allresp = zeros(length(data(plothese(1)).PID),length(plothese));
for i = 1:length(plothese)
	try
		allstim(:,i)=data(plothese(i)).PID(:)';
		allresp(:,i)=data(plothese(i)).ORN(:)';
	catch
		% pad with zeros
		allstim(:,i)=[0 data(plothese(i)).PID(:)'];
		allresp(:,i)=[0 data(plothese(i)).ORN(:)'];
		time = [NaN data(plothese(i)).time];
	end
end

plot(time,allstim)
set(gca,'XLim',[10 20])
xlabel('Time (s)')
ylabel('Odor concentration (V)')
title('Variability in stimulus presented in identical datasets')

PrettyFig;

snapnow;
delete(gcf);

%%
% The amplitude of the signal seems to vary significantly from dataset to dataset. However, the stimulus is well correlated: the following matrix shows the pairwise r-square between each piece of the data above:

r = NaN(length(plothese));
r2 = NaN(length(plothese));
for i = 1:length(r)-1
	for j = i+1:length(r)
		r(i,j) = rsquare(allstim(:,i),allstim(:,j));
		r2(i,j) = rsquare(allresp(:,i),allresp(:,j));
	end
end
disp(r)

%%
% How variable is the response of the neuron? The following plot shows the response to the neuron in these different datasets: 


figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(time,allresp)
set(gca,'XLim',[10 20])
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')
title('Variability in response presented in identical datasets')

PrettyFig;

snapnow;
delete(gcf);


%%
% The following matrix shows the pairwise r-square between each piece of the data above:

disp(r2)

%%
% So at this point it is unclear if the observed variability in stimulus is due to actual variability in stimulus amplitude,and the neuron adapts to this, or if the stimulus is actually identical, and the PID sensitivities are simply different from day to day. 

%%
% The following figure shows the results of backing the filters out of these supposedly identical datasets: 

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(filtertime,NormFilters(:,plothese))
xlabel('Lag (s)')
ylabel('Filter Amplitude (norm)')
set(gca,'XLim',[min(filtertime)-0.01 max(filtertime)+0.01])
set(gca,'YLim',[-0.9 1.2])
title('Variance in filter shape')

PrettyFig;

snapnow;
delete(gcf);

%%
% The following figure shows the results of the gain analysis on these supposedly identical datasets: 

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(history_lengths,low_slopes(:,plothese-1),'.-','Color',[0.5 1 0.5])
plot(history_lengths,high_slopes(:,plothese-1),'.-','Color',[1 0.5 0.5])
set(gca,'XScale','log')

xlabel('History Length (s)')
ylabel('Relative Gain')

% now plot the dots where significant
for i = (plothese-1)
	sig = p_values(:,i);
	sig = (sig<0.05);

	scatter(history_lengths(sig),low_slopes(sig,i),500,'g.')
	scatter(history_lengths(sig),high_slopes(sig,i),500,'r.')
end

PrettyFig;

snapnow;
delete(gcf);




return

%% 
% In the following plots, the response of the ORN is plotted against the best-fit LN model for an example history length. 

show_these = find(lowest_p < 0.05);
for i = 1:length(show_these)

	td = do_these(show_these(i));

	figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
	ph(3) = subplot(1,1,1);

	s = 300; % when we start for the gain analysis
	z = length(data(td).ORN) - 33; % where we end

	% compute the LNpred

	LinearFit = mean(data(td).ORN)+convolve(data(td).time,data(td).PID,Filters(:,show_these(i)),filtertime);
	LinearFit(LinearFit<0)=0;
	LNpred = hill(HillFit(:,show_these(i)),LinearFit);

	clear x
	x.response = data(td).ORN(s:z);
	x.prediction = LNpred(s:z);
	x.stimulus = data(td).PID(s:z);
	x.time = data(td).time(s:z);
	x.filter_length = 201;

	[~,loc]=min(p_values(:,show_these(i)));
	example_history_length = history_lengths(loc);

	GainAnalysis3(x,history_lengths,example_history_length,ph,NaN*history_lengths);

	title(strcat(data(td).neuron,'-',data(td).odor,'-\tau_H=',mat2str(example_history_length),'s' ))

	snapnow;
	delete(gcf);

end




