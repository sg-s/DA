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
N = length(data);
marker_size = 15;

% internal housekeeping: determine if being called by publish or not
calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end

if ~exist('HowGeneralIsGainAdaptation.mat','file')
	
	n = length(do_these);

	% initialise history lengths to run the analysis on.
	% history_lengths=[0:0.03:0.3 0.36:0.06:1 1.2:1.2:5];
	history_lengths=[0:3e-3:3e-2 0.036:3e-2:1 1.2:1.2:5];

	% initialise a matrix for all the linear filters
	Filters = NaN(200,N);

	% initialise a matrix for the parameters of all the non-linear functions
	HillFit = NaN(3,N);

	% initialise a matrix that stores all the slopes (and p-values) we calculate from the Gain Analysis
	low_slopes  = NaN(length(history_lengths),N);
	high_slopes = NaN(length(history_lengths),N);
	low_gof  = NaN(length(history_lengths),N);
	high_gof = NaN(length(history_lengths),N);
	p_values_low = NaN(length(history_lengths),N);
	p_values_high = NaN(length(history_lengths),N);
	all_slopes = NaN(length(history_lengths),N);
	data_min = NaN(length(history_lengths),N);
	data_max = NaN(length(history_lengths),N);
	low_min = NaN(length(history_lengths),N);
	low_max = NaN(length(history_lengths),N);
	high_min = NaN(length(history_lengths),N);
	high_max = NaN(length(history_lengths),N);

	% initialise a matrix that stores the r-square of the LN fit
	LNFitQuality = NaN(1,N);

	for i = 1:n

		td = do_these(i);
		disp(data(td).original_name)

		% fit a linear filter to data
		[K,~,filtertime] = FindBestFilter(data(td).PID,data(td).ORN,[],'filter_length=199;');
		Filters(:,td) = K;

		LinearFit = mean(data(td).ORN)+convolve(data(td).time,data(td).PID,K,filtertime);
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
		x = lsqcurvefit(@hill,[max(ydata) 2 2],xdata,ydata,[max(ydata)/2 2 1],[2*max(ydata) max(ydata) 10],fo);
		LNFit = hill(x,LinearFit);
		HillFit(:,td) = x;

		LNFitQuality(td) = rsquare(LNFit,data(td).ORN);
		if LNFitQuality(td) < rsquare(data(td).ORN,LinearFit)
			disp('WARNING: looks like the hill function fit failed')
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

		[p,low_slopes(:,td),high_slopes(:,td),low_gof(:,td),high_gof(:,td),~,extra_variables] = GainAnalysis4(x,history_lengths);
		p_values_low(:,td) = p(:,1);
		p_values_high(:,td) = p(:,2);
		data_min(:,td) = extra_variables.data_min(:);
		data_max(:,td) = extra_variables.data_max(:);
		low_min(:,td) = extra_variables.low_min(:);
		low_max(:,td) = extra_variables.low_max(:);
		high_min(:,td) = extra_variables.high_min(:);
		high_max(:,td) = extra_variables.high_max(:);
		all_slopes(:,td) = extra_variables.all_slopes(:);

		
	end

	% cache locally for use later
	filtertime = filtertime*mean(diff(data(td).time));
	save('HowGeneralIsGainAdaptation.mat','Filters','HillFit','LNFitQuality','high_slopes','low_slopes','p_values_low','p_values_high','filtertime','history_lengths','low_gof','high_gof','all_slopes','data_min','data_max','low_min','low_max','high_min','high_max')
else
	load('HowGeneralIsGainAdaptation.mat')
end


% corrects for points where the p-value is reported as 0 by the bootstrap. it can't be, it's just too small to measure. 
p_values_low(p_values_low==0) = min(nonzeros(p_values_low(:)));
p_values_high(p_values_high==0) = min(nonzeros(p_values_high(:)));

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

if being_published
	snapnow;
	delete(gcf);
end


%       ##     ## #### ##       ##          ######## ##     ## ##    ##  ######  
%       ##     ##  ##  ##       ##          ##       ##     ## ###   ## ##    ## 
%       ##     ##  ##  ##       ##          ##       ##     ## ####  ## ##       
%       #########  ##  ##       ##          ######   ##     ## ## ## ## ##       
%       ##     ##  ##  ##       ##          ##       ##     ## ##  #### ##       
%       ##     ##  ##  ##       ##          ##       ##     ## ##   ### ##    ## 
%       ##     ## #### ######## ########    ##        #######  ##    ##  ######  

%% Variance of Output Non-linearities 
% In each data set, after backing out the filters, we then fit a 3-parameter hill function to the residuals. The following function shows the shape of all the Hill functions. Each curve is the best fit for one data set. 

y = NaN(251,N);
x = NaN(251,N);

for i = do_these
	xx = 0:1:max(data(i).ORN);
	xx = [NaN(1,length(y)-length(xx)) xx];
	y(:,i) = hill(HillFit(:,i),xx);
	x(:,i) = xx;
end

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(x,y)
xlabel('Linear Output (Hz)')
ylabel('Nonlinearity Output (Hz)')
set(gca,'XLim',[0 255])
title('Variance in Nonlinearity shape')

PrettyFig;


if being_published
	snapnow;
	delete(gcf);
end


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
for i = do_these
	sig_low = p_values_low(:,i);
	sig_high = p_values_high(:,i);
	sig_low = (sig_low<0.05);
	sig_high = (sig_high<0.05);

	scatter(history_lengths(sig_low),low_slopes(sig_low,i),500,'g.')
	scatter(history_lengths(sig_high),high_slopes(sig_high,i),500,'r.')
end

xlabel('History Length (s)')
ylabel('Relative Gain')
PrettyFig;

set(gca,'XLim',[min(nonzeros(history_lengths))/2 2*max(history_lengths)])
if being_published
	snapnow;
	delete(gcf);
end


%%
% However, when we go to very short history lengths, we end up sampling only points in time when the neuron is not firing at all. These are pathological points and cause wildly inaccurate estimates of gain for the low stimuli. To see what we are talking about, the following figure shows the maximum of the subset of the data sampled, normalised by the maximum of the data, as a function of history length for the various data sets analysed here. 

figure('outerposition',[0 0 700 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(history_lengths,low_max./data_max)
set(gca,'XScale','log','YLim',[-0.1 1.1])
xlabel('History Length (s)')
ylabel('Max of "Low" Response Data')

PrettyFig;

if being_published
	snapnow;
	delete(gcf);
end

%%
% We arbitrarily decide to throw out all data where the sampled subset is less than half the full data set. 

low_slopes((low_max./data_max)<0.5) = NaN;
p_values_low((low_max./data_max)<0.5) = NaN;


% %%
% % However, in some cases, the fits to the clouds of points in the gain analysis may not be very good. The following plot shows the distribution of r-square values of the fits that are used to determine the gain in each of these cases, for the entire dataset:

% figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
% [y,x]=hist(low_gof(:),30);
% plot(x,y,'g')
% [y,x]=hist(high_gof(:),30);
% plot(x,y,'r')
% legend({'Low Slopes','High Slopes'})
% xlabel('r-square of fit')
% ylabel('Count')
% title('Some fits are very poor')

% PrettyFig;

% if being_published
% 	snapnow;
% 	delete(gcf);
% end


% %%
% % If we retain only the points where the r-square of the fit is >0.8, we end up retaining the following percent of the low slopes:

% disp(100*length(low_gof(low_gof>0.8))/(length(low_gof(:))))

% %%
% % and the following percent of the high-slopes data:

% disp(100*length(high_gof(high_gof>0.8))/(length(high_gof(:))))

%%
% In the following plot, we only retain this data:

% low_slopes(low_gof<0.8) = NaN;
% high_slopes(high_gof<0.8) = NaN;

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(history_lengths,low_slopes,'.-','Color',[0.5 1 0.5])
plot(history_lengths,high_slopes,'.-','Color',[1 0.5 0.5])
set(gca,'XScale','log')

xlabel('History Length (s)')
ylabel('Relative Gain')

% now plot the dots where significant
for i = do_these
	sig_low = p_values_low(:,i);
	sig_high = p_values_high(:,i);
	sig_low = (sig_low<0.05);
	sig_high = (sig_high<0.05);

	scatter(history_lengths(sig_low),low_slopes(sig_low,i),500,'g.')
	scatter(history_lengths(sig_high),high_slopes(sig_high,i),500,'r.')
end

set(gca,'XLim',[1e-3 6])

PrettyFig;

if being_published
	snapnow;
	delete(gcf);
end


%%
% In general, are the green slopes (gain following low stimuli) significantly higher than the red slopes (gain following high stimuli)? To determine this, we plot the slope of the best fit line vs. the p value of that fit. 

figure('outerposition',[0 0 1100 600],'PaperUnits','points','PaperSize',[1100 600]); hold on
subplot(1,2,1), hold on
d = low_slopes(:)./all_slopes(:); %  - high_slopes(:);
plot(p_values_low(:),d,'g.','MarkerSize',marker_size) 
set(gca,'XScale','log')
% plot a reference line for equal gain
plot([.01 1],[1 1 ],'k--')
% plot a reference line for p = 0.05
plot([.05 .05],[min(d)/2 max(d)*2 ],'k--')
set(gca,'XLim',[min(p_values_low(:))/2 1.1],'YLim',[min(d)/2 max(d)*1.1])
ylabel('Relative Gain (low stimuli)')
xlabel('p-value')

subplot(1,2,2), hold on
d = high_slopes(:)./all_slopes(:); %  - high_slopes(:);
plot(p_values_high(:),d,'r.','MarkerSize',marker_size) 
set(gca,'XScale','log')
% plot a reference line for equal gain
plot([.01 1],[1 1 ],'k--')
% plot a reference line for p = 0.05
plot([.05 .05],[min(d)/2 max(d)*2 ],'k--')
set(gca,'XLim',[min(p_values_high(:))/2 1.1],'YLim',[min(d)/2 max(d)*1.1])
ylabel('Relative Gain (high stimuli)')
xlabel('p-value')

PrettyFig;

if being_published
	snapnow;
	delete(gcf);
end



%% 
% The reference horizontal line indicates equal slope (equal gain for low and high stimuli) and the reference vertical line indicates a p-value of _p=0.05_. Of the statistically significant differences in gain (left half-plane), there are far more points in the top-left quadrant (gain enhancement to low stimuli) than there are in the bottom-left quadrant (gain suppression to low stimuli). 



%%
% Every data set tested showed statistically significant gain control at some timescale. 

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
allstim = zeros(length(data(plothese(1)).PID),N);
allresp = zeros(length(data(plothese(1)).PID),N);
for i = plothese
	try
		allstim(:,i)=data(i).PID(:)';
		allresp(:,i)=data(i).ORN(:)';
	catch
		% pad with zeros
		allstim(:,i)=[data(i).PID(:)' 0];
		allresp(:,i)=[data(i).ORN(:)' 0];
		time = [NaN data(i).time];
	end
end

plot(time,allstim)
set(gca,'XLim',[10 20])
xlabel('Time (s)')
ylabel('Odor concentration (V)')
title('Variability in stimulus presented in identical datasets')

PrettyFig;

if being_published
	snapnow;
	delete(gcf);
end


%%
% The amplitude of the signal seems to vary significantly from dataset to dataset. However, the stimulus is well correlated: the following matrix shows the pairwise r-square between each piece of the data above:

r = NaN(length(plothese));
r2 = NaN(length(plothese));
for i = 1:length(r)-1
	a = plothese(i);
	for j = i+1:length(r)
		b = plothese(j);
		r(i,j) = rsquare(allstim(:,a),allstim(:,b));
		r2(i,j) = rsquare(allresp(:,a),allresp(:,b));
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

if being_published
	snapnow;
	delete(gcf);
end


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

if being_published
	snapnow;
	delete(gcf);
end

%%
% The following figure shows the results of the gain analysis on these supposedly identical datasets: 

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(history_lengths,low_slopes(:,plothese),'.-','Color',[0.5 1 0.5])
plot(history_lengths,high_slopes(:,plothese),'.-','Color',[1 0.5 0.5])
set(gca,'XScale','log')

xlabel('History Length (s)')
ylabel('Relative Gain')

% now plot the dots where significant
for i = plothese
	sig_low = p_values_low(:,i);
	sig_low = (sig_low<0.05);
	sig_high = p_values_high(:,i);
	sig_high = (sig_high<0.05);

	scatter(history_lengths(sig_low),low_slopes(sig_low,i),500,'g.')
	scatter(history_lengths(sig_high),high_slopes(sig_high,i),500,'r.')
end

set(gca,'XLim',[min(nonzeros(history_lengths))/2 max(history_lengths)*2])
PrettyFig;

if being_published
	snapnow;
	delete(gcf);
end



%  ######   #######  ########  ########     ##       ######## ##    ##  ######   ######## ##     ## 
% ##    ## ##     ## ##     ## ##     ##    ##       ##       ###   ## ##    ##     ##    ##     ## 
% ##       ##     ## ##     ## ##     ##    ##       ##       ####  ## ##           ##    ##     ## 
% ##       ##     ## ########  ########     ##       ######   ## ## ## ##   ####    ##    ######### 
% ##       ##     ## ##   ##   ##   ##      ##       ##       ##  #### ##    ##     ##    ##     ## 
% ##    ## ##     ## ##    ##  ##    ##     ##       ##       ##   ### ##    ##     ##    ##     ## 
%  ######   #######  ##     ## ##     ##    ######## ######## ##    ##  ######      ##    ##     ## 

%% Effect of correlation length
% What is the effect of the correlation length of the stimulus on the estimation of the degree of gain control? The following subset of the data is analysed, where the same odor is presented to the same type of neuron, but with different correlation lengths. All the data was recorded on the same day:


plothese = [3     4   5];
tau_c	 = [100   30 50];
tau_c_text = {'100ms','30ms','50ms'};
for i = 1:length(plothese)
	disp(data(plothese(i)).original_name)
end

%%
% The following figure shows that the autocorrelation times of the stimulus are indeed different:

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
t = 0:3e-3:3e-1;
a = NaN(length(t),length(plothese));
ti=1;
for i = plothese
	a(:,ti)=autocorr(data(i).PID,100);
	ti = ti+1;
end
plot(t,a)
legend(tau_c_text)
xlabel('Lag (s)')
ylabel('Autocorrelation')

PrettyFig;

if being_published
	snapnow;
	delete(gcf);
end



%%
% The following figure shows the results of the gain analysis on these datasets: 

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(history_lengths,low_slopes(:,plothese(1)),'.-','Color',[0 1.0 0])
plot(history_lengths,low_slopes(:,plothese(2)),'.-','Color',[0 0.4 0])
plot(history_lengths,low_slopes(:,plothese(3)),'.-','Color',[0 0.7 0])

plot(history_lengths,high_slopes(:,plothese(1)),'.-','Color',[1.0 0 0])
plot(history_lengths,high_slopes(:,plothese(2)),'.-','Color',[0.4 0 0])
plot(history_lengths,high_slopes(:,plothese(3)),'.-','Color',[0.7 0 0])
set(gca,'XScale','log')

xlabel('History Length (s)')
ylabel('Relative Gain')

% now plot the dots where significant
% now plot the dots where significant
for i = plothese
	sig_low = p_values_low(:,i);
	sig_low = (sig_low<0.05);
	sig_high = p_values_high(:,i);
	sig_high = (sig_high<0.05);

	scatter(history_lengths(sig_low),low_slopes(sig_low,i),500,'g.')
	scatter(history_lengths(sig_high),high_slopes(sig_high,i),500,'r.')
end

set(gca,'XLim',[min(nonzeros(history_lengths))/2 max(history_lengths)*2])
legend(tau_c_text)

PrettyFig;

if being_published
	snapnow;
	delete(gcf);
end


%%
% So it looks like the peak of the green curves (gain in response to low stimuli) grows smaller, the smaller the correlation length of the stimulus. 

%   ########  #### ######## ########    ##    ## ######## ##     ## ########   #######  ##    ## 
%   ##     ##  ##  ##       ##          ###   ## ##       ##     ## ##     ## ##     ## ###   ## 
%   ##     ##  ##  ##       ##          ####  ## ##       ##     ## ##     ## ##     ## ####  ## 
%   ##     ##  ##  ######   ######      ## ## ## ######   ##     ## ########  ##     ## ## ## ## 
%   ##     ##  ##  ##       ##          ##  #### ##       ##     ## ##   ##   ##     ## ##  #### 
%   ##     ##  ##  ##       ##          ##   ### ##       ##     ## ##    ##  ##     ## ##   ### 
%   ########  #### ##       ##          ##    ## ########  #######  ##     ##  #######  ##    ## 

%% Effect of receptor/ORN
% In the following section, we analyse a subset of the data where the same odor is presented in the same manner, but to two different neurons: ab3A and pb1A. The data used is:

plothese = [18 12];
neuron	 = {'ab3A','pb1A'};
for i = 1:length(plothese)
	disp(data(plothese(i)).original_name)
end



%%
% The following figure shows the two stimulus in these two cases:

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(data(plothese(1)).time,data(plothese(1)).PID,'k')
plot(data(plothese(2)).time,data(plothese(2)).PID,'r')
set(gca,'XLim',[10 20])
xlabel('Time (s)')
ylabel('Odor concentration (V)')
legend(neuron);
PrettyFig;

if being_published
	snapnow;
	delete(gcf);
end

%%
% It looks like one stimulus is massively bigger than the other. Normalising by the mean, we get:
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(data(plothese(1)).time,data(plothese(1)).PID/mean(data(plothese(1)).PID),'k')
plot(data(plothese(2)).time,data(plothese(2)).PID/mean(data(plothese(2)).PID),'r')
set(gca,'XLim',[10 20])
xlabel('Time (s)')
ylabel('Odor concentration (norm)')
legend(neuron);
PrettyFig;

if being_published
	snapnow;
	delete(gcf);
end

%%
% The responses of the two different neurons to this temporally identical stimulus are:

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(data(plothese(1)).time,data(plothese(1)).ORN,'k')
plot(data(plothese(2)).time,data(plothese(2)).ORN,'r')
set(gca,'XLim',[10 20])
xlabel('Time (s)')
ylabel('Neuron Response (Hz)')
legend(neuron);

PrettyFig;

snapnow;
delete(gcf);


%%
% The following figure shows how the gain analysis of these two datasets differ:

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(history_lengths,low_slopes(:,plothese(1)),'.-','Color',[0 1.0 0])
plot(history_lengths,low_slopes(:,plothese(2)),'.-','Color',[0 0.4 0])

plot(history_lengths,high_slopes(:,plothese(1)),'.-','Color',[1.0 0 0])
plot(history_lengths,high_slopes(:,plothese(2)),'.-','Color',[0.4 0 0])

set(gca,'XScale','log')

xlabel('History Length (s)')
ylabel('Relative Gain')

% now plot the dots where significant
for i = plothese
	sig_low = p_values_low(:,i);
	sig_low = (sig_low<0.05);
	sig_high = p_values_high(:,i);
	sig_high = (sig_high<0.05);

	scatter(history_lengths(sig_low),low_slopes(sig_low,i),500,'g.')
	scatter(history_lengths(sig_high),high_slopes(sig_high,i),500,'r.')
end

legend(neuron)
set(gca,'XLim',[min(nonzeros(history_lengths))/2 max(history_lengths)*2])
PrettyFig;

if being_published
	snapnow;
	delete(gcf);
end



%       ########  #### ######## ########         #######  ########   #######  ########  
%       ##     ##  ##  ##       ##              ##     ## ##     ## ##     ## ##     ## 
%       ##     ##  ##  ##       ##              ##     ## ##     ## ##     ## ##     ## 
%       ##     ##  ##  ######   ######          ##     ## ##     ## ##     ## ########  
%       ##     ##  ##  ##       ##              ##     ## ##     ## ##     ## ##   ##   
%       ##     ##  ##  ##       ##       ###    ##     ## ##     ## ##     ## ##    ##  
%       ########  #### ##       ##       ###     #######  ########   #######  ##     ## 

% %% Effect of odor 
% In this section, we investigate the effect of varying the odor type, while keeping the temporal structure and the type of neuron it is presented to fixed. We use the following datasets, where the valve is switched with the same temporal pattern, and the stimulus is presented to ab3A neurons:

% plothese = [2 4 7];
% odor = {data(plothese).odor};

% for i = 1:length(plothese)
% 	disp(data(plothese(i)).original_name)
% end

% %%
% % The following figure shows the stimulus from these datasets, normalised by the mean:

% figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
% plot(data(plothese(1)).time,data(plothese(1)).PID/mean(data(plothese(1)).PID),'b')
% plot(data(plothese(2)).time,data(plothese(2)).PID/mean(data(plothese(2)).PID),'g')
% plot(data(plothese(3)).time,data(plothese(3)).PID/mean(data(plothese(3)).PID),'r')
% set(gca,'XLim',[10 20])
% xlabel('Time (s)')
% ylabel('Odor concentration (norm)')
% legend(odor);
% PrettyFig;

% if being_published
% 	snapnow;
% 	delete(gcf);
% end

% %%
% % The responses of ab3A to these three different odors looks like. The panel on the left shows the actual responses, and the panel on the right shows the normalised histogram of the responses of the neurons.

% figure('outerposition',[0 0 1200 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
% subplot(1,5,1:4), hold on
% plot(data(plothese(1)).time,data(plothese(1)).ORN,'Color','b')
% plot(data(plothese(2)).time,data(plothese(2)).ORN,'Color','g')
% plot(data(plothese(3)).time,data(plothese(3)).ORN,'Color','r')
% set(gca,'XLim',[10 20])
% xlabel('Time (s)')
% ylabel('Firing rate (Hz)')
% legend(odor);

% subplot(1,5,5), hold on
% [x,y] = hist(data(plothese(1)).ORN,30); plot(x/max(x),y,'b')
% [x,y] = hist(data(plothese(2)).ORN,30); plot(x/max(x),y,'g')
% [x,y] = hist(data(plothese(3)).ORN,30); plot(x/max(x),y,'r')

% PrettyFig;

% if being_published
% 	snapnow;
% 	delete(gcf);
% end

% %%
% % We now look at how the gain analysis of these datasets differs: 

% figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
% plot(history_lengths,low_slopes(:,plothese(1)),'.-','Color',[0 1.0 0])
% plot(history_lengths,low_slopes(:,plothese(2)),'.-','Color',[0 0.7 0])
% plot(history_lengths,low_slopes(:,plothese(3)),'.-','Color',[0 0.4 0])

% plot(history_lengths,high_slopes(:,plothese(1)),'.-','Color',[1.0 0 0])
% plot(history_lengths,high_slopes(:,plothese(2)),'.-','Color',[0.7 0 0])
% plot(history_lengths,high_slopes(:,plothese(3)),'.-','Color',[0.4 0 0])

% set(gca,'XScale','log','XLim',[1e-3 5])

% xlabel('History Length (s)')
% ylabel('Relative Gain')

% % now plot the dots where significant
% for i = plothese
% 	sig_low = p_values_low(:,i);
% 	sig_low = (sig_low<0.05);
% 	sig_high = p_values_high(:,i);
% 	sig_high = (sig_high<0.05);

% 	scatter(history_lengths(sig_low),low_slopes(sig_low,i),500,'g.')
% 	scatter(history_lengths(sig_high),high_slopes(sig_high,i),500,'r.')
% end

% legend(odor)

% PrettyFig;

% if being_published
% 	snapnow;
% 	delete(gcf);
% end


%% Gain Analysis Examples
% In the following plots, the response of the ORN is plotted against the best-fit LN model for an example history length. 

for i = do_these



	s = 300; % when we start for the gain analysis
	z = length(data(i).ORN) - 33; % where we end

	% compute the LNpred
	LinearFit = mean(data(i).ORN)+convolve(data(i).time,data(i).PID,Filters(:,i),filtertime);
	LinearFit(LinearFit<0)=0;
	LNpred = hill(HillFit(:,i),LinearFit);

	clear x
	x.response = data(i).ORN(s:z);
	x.prediction = LNpred(s:z);
	x.stimulus = data(i).PID(s:z);
	x.time = data(i).time(s:z);
	x.filter_length = 201;


	r = low_slopes(:,i) - high_slopes(:,i);
	r(p_values_low(:,i)>0.05) = 0;
	r(p_values_high(:,i)>0.05) = 0;
	[~,loc]=max(r);

	if max(r) > 0
		% ok, we can show this
		figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
		ph(3) = subplot(1,1,1);

		example_history_length = history_lengths(loc);

		GainAnalysis3(x,history_lengths,example_history_length,ph,NaN*history_lengths);

		title(strcat(data(i).neuron,'-',data(i).odor,'-\tau_H=',mat2str(example_history_length),'s' ))

		if being_published
			snapnow;
			delete(gcf);
		end

	else
		% can't show this

	end


end


%% Version Info
% The file that generated this document is called:
disp(mfilename)

%%
% and its md5 hash is:
Opt.Input = 'file';
disp(DataHash(strcat(mfilename,'.m'),Opt))



