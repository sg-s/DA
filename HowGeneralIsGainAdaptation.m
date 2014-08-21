% HowGeneralIsGainAdaptation.m
% summarises gain adaptation in all of carlotta's data. makes plots organised by receptor, odor, correlation length, etc. 
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


% the following section pre-computes the processor-heavy parts and caches them for later use by publish()
load('/local-data/DA-paper/data.mat')

if ~exist('HowGeneralIsGainAdaptation.mat','file')
	do_these = 2:21;
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
	p_values = NaN(length(history_lengths),n);

	% initialise a matrix that stores the r-square of the LN fit
	LNFitQuality = NaN(1,n);

	for i = 1:n

		td = do_these(i);

		% fit a linear filter to data
		[K,~,filtertime] = FindBestFilter(data(td).PID,data(td).ORN,[],'filter_length=199;');
		Filters(:,i) = K;

		LinearFit = mean(data(td).ORN)+convolve(data(td).time,data(td).PID,data(td).K,data(td).filtertime);
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
		if LNFitQuality(i) < 0.9
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

		[p_values(:,i),low_slopes(:,i),high_slopes(:,i)] = GainAnalysis3(x,history_lengths);
		
	end

	% cache locally for use later
	filtertime = filtertime*mean(diff(data(td).time));
	save('HowGeneralIsGainAdaptation.mat','Filters','HillFit','LNFitQuality','high_slopes','low_slopes','p_values','filtertime')
else
	load('HowGeneralIsGainAdaptation.mat')
end


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

PrettyFig;

snapnow;
delete(gcf);

%% Variance of Output Nonlinearities 
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


PrettyFig;

%% Summary of Gain Analysis in all data

snapnow;
delete(gcf);