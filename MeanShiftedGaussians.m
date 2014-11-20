% 
% 
% created by Srinivas Gorur-Shandilya at 8:38 , 12 November 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


% internal housekeeping: determine if being called by publish or not
calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end
% redo = 1; deliberately unset

data_root = '/local-data/DA-paper/fast-flicker/orn/';
allfiles  = dir(strcat(data_root,'*.mat'));

% combine all data
% if redo
% 	combined_data = ReduceORNData(data_root,allfiles);
% end

% paradigm_names = unique(combined_data.paradigm);

% load cached data
load('MeanShiftedGaussians.mat')


%% Mean Shifted Gaussians
% How to ORNs respond to mean shifted gaussians? Specifcally, how do they vary their input-output curve? Is the adaptation to this mean optimal (a la Laughlin etc)? Can we find evidence for fast gain adaptation in the curves themselves? 

%         ######  ######## #### ##     ## ##     ## ##       ##     ##  ######  
%        ##    ##    ##     ##  ###   ### ##     ## ##       ##     ## ##    ## 
%        ##          ##     ##  #### #### ##     ## ##       ##     ## ##       
%         ######     ##     ##  ## ### ## ##     ## ##       ##     ##  ######  
%              ##    ##     ##  ##     ## ##     ## ##       ##     ##       ## 
%        ##    ##    ##     ##  ##     ## ##     ## ##       ##     ## ##    ## 
%         ######     ##    #### ##     ##  #######  ########  #######   ######  

%% Stimulus Distributions 
% Three gaussians with different means are chosen. The following figure shows the distributions for every trial of the data analysed here, colour coded by the experimental paradigm. 

% some global parameters
nbins = 50;
histx = [];
histy = [];
paradigm = [];
dt = 3e-3;
all_pid = [];

a = floor(30/3e-3);
z = floor(55/3e-3);

% assemble all histograms
for i = 1:length(paradigm_names)
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	for j = 1:length(plot_these)
		this_pid = combined_data.PID(plot_these(j),:);
		this_pid = this_pid(a:z);
		[y,x] = hist(this_pid,nbins);
		if mean(x) > .5
			histx = [histx; x];
			histy = [histy; y];
			paradigm = [paradigm i];
			all_pid = [all_pid; this_pid];
		end
	end
end

					



c = parula(length(paradigm_names));
paradigm = paradigm -  min(paradigm);
paradigm = paradigm  + 1;

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
for i = 1:length(paradigm)
	plot(histx(i,:),histy(i,:),'Color',c(paradigm(i),:))
end

xlabel('PID (V)')
ylabel('Count')
titlestr = strcat(mat2str(length(paradigm)),' trials');
title(titlestr)


PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%% Stimulus Reproducibility 
% In this section, we look at the reproducibility of the stimulus. The following figure shows the stimulus for all the data we look at here, plotted on top of each other, colour-coded by experimental paradigm. 

c = parula(length(paradigm_names));
time = 1:length(all_pid);
time = time*dt;

figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
for i = 1:width(all_pid)
	plot(time,all_pid(i,:),'Color',c(paradigm(i),:))
end

ylabel('PID (V)')
xlabel('Time (s)')
titlestr = strcat(mat2str(length(paradigm)),' trials');
title(titlestr)

set(gca,'XLim',[15 20])

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%          ########  ########  ######  ########   #######  ##    ##  ######  ######## 
%          ##     ## ##       ##    ## ##     ## ##     ## ###   ## ##    ## ##       
%          ##     ## ##       ##       ##     ## ##     ## ####  ## ##       ##       
%          ########  ######    ######  ########  ##     ## ## ## ##  ######  ######   
%          ##   ##   ##             ## ##        ##     ## ##  ####       ## ##       
%          ##    ##  ##       ##    ## ##        ##     ## ##   ### ##    ## ##       
%          ##     ## ########  ######  ##         #######  ##    ##  ######  ######## 

%% Neuron Responses: Overview
% The following figure shows the responses of the ORNs to this stimuli, and their distributions. 

figure('outerposition',[0 0 1400 500],'PaperUnits','points','PaperSize',[1400 500]); hold on
subplot(1,5,1:4), hold on
for i = 1:length(paradigm_names)
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	this_orn=mean2(combined_data.fA(:,plot_these));
	time = 3e-3*(1:length(this_orn));
	plot(time,this_orn)
end
ylim = get(gca,'YLim');
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')

subplot(1,5,5), hold on
for i = 1:length(paradigm_names)
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	this_orn=mean2(combined_data.fA(:,plot_these));
	this_orn = this_orn(a:z);
	[x,y] = hist(this_orn,50);
	plot(x,y)

end
set(gca,'YLim',ylim);
xlabel('Count')
title('Distribution')


PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%% Correlation Times 
% On what timescales are the stimulus and the response correlated? In the following figure, we plot the autocorrelation function of the stimulus and the responses, for each stimulus presented:

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on

c = parula(length(paradigm_names));

subplot(1,2,1), hold on
for i = 1:length(paradigm_names)
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	this_pid=mean2(combined_data.PID(plot_these,:));
	this_pid= this_pid(a:z);
	[y,x]=autocorr(this_pid,500);
	x = x*3e-3;
  	plot(x,y,'Color',c(i,:))
end
set(gca,'XScale','log','XMinorTick','on')
xlabel('Time (s)')
ylabel('Autocorrelation')
title('Stimulus')


subplot(1,2,2), hold on
for i = 1:length(paradigm_names)
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	this_orn=mean2(combined_data.fA(:,plot_these));
	this_orn = this_orn(a:z);
	[y,x]=autocorr(this_orn,500);
	x = x*3e-3;
  	plot(x,y,'Color',c(i,:))
end
set(gca,'XScale','log','XMinorTick','on')
xlabel('Time (s)')
ylabel('Autocorrelation')
title('ORN Responses')

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%% Trends in Data
% Are there any trends in the data? In the following figure, we coarse-grain the data by binning everything along 5-second bins to look at long-term trends in the data. The various colors correspond to various stimulus means, and correspond to other figures in this document. 

% plot_data is indexed by where we start
all_start = [15:5:50];
all_end = all_start+5;


for i = 1:length(paradigm_names)
	plot_data(i).stim_slope = [];
	plot_data(i).stim_slope_err = [];
	plot_data(i).stim_mean = [];
	plot_data(i).stim_mean_err = [];
	plot_data(i).resp_slope = [];
	plot_data(i).resp_slope_err = [];
	plot_data(i).resp_mean = [];
	plot_data(i).resp_mean_err = [];

	for j = 1:length(all_start)
		a = floor(all_start(j)/3e-3);
		z = floor(all_end(j)/3e-3);
		n = sqrt(z-a);

		plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
		these_pid=mean2(combined_data.PID(plot_these,:));
		these_resp=mean2(combined_data.fA(:,plot_these));

		cropped_pid = these_pid(a:z);
		cropped_resp = these_resp(a:z);

		plot_data(i).stim_mean = 		[plot_data(i).stim_mean mean(cropped_pid)];
		plot_data(i).stim_mean_err = 	[plot_data(i).stim_mean_err std(cropped_pid)/n];

		plot_data(i).resp_mean = 		[plot_data(i).resp_mean mean(cropped_resp)];
		plot_data(i).resp_mean_err = 	[plot_data(i).resp_mean_err std(cropped_resp)/n];


	end
end

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
c = parula(length(paradigm_names));
for i = 1:length(plot_data)
	errorbar(all_start+2.5,plot_data(i).stim_mean,plot_data(i).stim_mean_err)

end
xlabel('Time (s)')
ylabel('PID (V)')

subplot(1,2,2), hold on
c = parula(length(paradigm_names));
for i = 1:length(plot_data)
	errorbar(all_start+2.5,plot_data(i).resp_mean,plot_data(i).resp_mean_err)

end
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')


PrettyFig;
if being_published
	snapnow
	delete(gcf)
end


return


%% Neuron Responses: Input-Output Curve Changes
% How do ORNs change their input-output curves when presented with these different stimuli? If the ORN is a perfect sensor, that adapts perfectly, it would move and stretch its I/O curve so that it is equal to cumulative density function of the stimulus. However, we know that ORNs don't adapt perfectly, so they must do something else. What is it? If the change dominated by a lateral movement of a stretch?

% this stores the LN model parameters for all the data
if redo
	clear LNModel
	LNModel.K = [];
	LNModel.H = []; % stores hill parameters
	LNModel.LinearFit = [];
	LNModel.LNFit = [];
	LNModel.H_domain = [];


	for i = 1:length(paradigm_names)
		plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
		this_orn=mean2(combined_data.fA(:,plot_these));
		this_pid = mean2(combined_data.PID(plot_these,:));

		time = 3e-3*(1:length(this_orn));

		this_pid = this_pid(a:z);
		this_orn = this_orn(a:z);
		time = time(a:z);

		% fit a polynomial trend
		ff = fit(time(:),this_pid(:),'poly2');
		this_pid = this_pid - ff(time)';

		% normalise stimulus and response
		% this_pid = this_pid - mean(this_pid);
		% this_pid = this_pid/std(this_pid);
		% this_orn = this_orn - mean(this_orn);
		% this_orn = this_orn/std(this_orn);

		[K,~,filtertime] = FindBestFilter(this_pid,this_orn,[],'filter_length=299;');
		filtertime = filtertime*mean(diff(time));
		K = K/max(K);
		LinearFit = convolve(time,this_pid,K,filtertime) + mean(this_orn);

		xdata = LinearFit(:);
		ydata = this_orn(:);

		ydata(isnan(xdata)) = [];
		xdata(isnan(xdata)) = [];


		fo=optimset('MaxFunEvals',1000,'Display','none');
		x = lsqcurvefit(@hill4,[max(ydata) 2 2 0],xdata,ydata,[max(ydata)/100 -max(ydata) 1 0],[20*max(ydata) 10*max(ydata) 10 max(ydata)],fo);
		LNFit = hill4(x,xdata);

		% save all of this for later
		LNModel(i).K = K;
		LNModel(i).H = x;
		LNModel(i).LinearFit = LinearFit;
		LNModel(i).LNFit = LNFit;
		LNModel(i).H_domain = xdata;
		LNModel(i).H_range =  ydata;
		LNModel(i).this_orn = this_orn;
		LNModel(i).LNFit_r2 = rsquare(LNFit,ydata);
		LNModel(i).LinearFit_r2 = rsquare(xdata,ydata);
	end

end


%%
% The following figure shows each of the ORN responses to each stimulus together with the best linear fits:


figure('outerposition',[0 0 1400 600],'PaperUnits','points','PaperSize',[1400 600]); hold on
for i = 1:6
	subplot(2,3,i), hold on
	this_orn = LNModel(i).this_orn;
	time = 3e-3*(1:length(this_orn));
	time = time+a*3e-3;
	plot(time,this_orn,'k');
	time = 3e-3*(1:length(LNModel(i).LinearFit));
	time = time+a*3e-3;
	lh=plot(time,LNModel(i).LinearFit,'r');
	set(gca,'XLim',[30 40])
	legend(lh,oval(LNModel(i).LinearFit_r2,2))
	title(paradigm_names{i})
end 
PrettyFig;


if being_published
	snapnow
	delete(gcf)
end

%%
% The following figure shows each of the ORN responses to each stimulus together with the best LN Model fits:


figure('outerposition',[0 0 1400 600],'PaperUnits','points','PaperSize',[1400 600]); hold on
for i = 1:6
	subplot(2,3,i), hold on
	this_orn = LNModel(i).this_orn;
	time = 3e-3*(1:length(this_orn));
	time = time+a*3e-3;
	plot(time,this_orn,'k');
	time = 3e-3*(1:length(LNModel(i).LNFit));
	time = time+a*3e-3;
	lh=plot(time,LNModel(i).LNFit,'r');
	set(gca,'XLim',[30 40])
	legend(lh,oval(LNModel(i).LNFit_r2,2))
	title(paradigm_names{i})
end 
PrettyFig;


if being_published
	snapnow
	delete(gcf)
end


%%
% The following figure shows the filters backed out of each data set (corresponding to each mean), and the corresponding non-linearities.  

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
c = parula(length(LNModel));
for i = 1:length(LNModel)
	plot(filtertime,LNModel(i).K,'Color',c(i,:))
end
xlabel('Filter Lag (s)')
ylabel('Filter Amplitude (norm)')

subplot(1,2,2), hold on
for i = 1:length(LNModel)
	x = sort(LNModel(i).H_domain);
	plot(x,hill4(LNModel(i).H,x),'Color',c(i,:))
end
xlabel('Filter Output (a.u.)')
ylabel('ORN Firing Rate (Hz)')

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end


%%
% Are these changes largely subtractive or largely divisive? In the following figure, we attempt to fit the data using only subtractive, only divisive, or mixed models. 

lh= [];
L = {};
figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1), hold on
for i = 1:length(LNModel)
	x = sort(LNModel(i).H_domain);
	lh = [lh plot(x,hill4(LNModel(i).H,x),'Color',c(i,:))];
	L = [L oval(LNModel(i).LNFit_r2,2)];
end
title('Mixed Model')
xlabel('Filtered Stimulus (a.u.)')
ylabel('Static Nonlinearity (Hz)')
legend(lh,L,'Location','southeast')

% fixed parameters
Kd0 = LNModel(1).H(2);
n0 = LNModel(1).H(3);

subplot(1,3,2), hold on
lh= [];
L = {};
title('Divisive Only')
for i = 1:length(LNModel)
	xdata = LNModel(i).H_domain;
	ydata = LNModel(i).H_range;

	fo=optimset('MaxFunEvals',1000,'Display','none');
	x = lsqcurvefit(@hill4,[max(ydata) Kd0 2 0],xdata,ydata,[max(ydata)/2 Kd0-1e-6 1 0],[2*max(ydata) Kd0+1e-6 10 max(ydata)],fo);
	LNFit = hill4(x,xdata);

	xdata = sort(xdata);

	lh = [lh plot(xdata,hill4(x,xdata),'Color',c(i,:))];
	r2 = rsquare(LNFit,ydata);
	L = [L oval(r2,2)];
end
xlabel('Filtered Stimulus (a.u.)')
ylabel('Static Nonlinearity (Hz)')
legend(lh,L,'Location','southeast')


subplot(1,3,3), hold on
lh= [];
L = {};
title('Subtractive Only')
for i = 1:length(LNModel)
	xdata = LNModel(i).H_domain;
	ydata = LNModel(i).H_range;

	fo=optimset('MaxFunEvals',1000,'Display','none');
	x = lsqcurvefit(@hill4,[max(ydata) 10 n0 0],xdata,ydata,[max(ydata)/2 0 n0-1e-6 0],[2*max(ydata) max(ydata) n0+1e-6 max(ydata)],fo);
	LNFit = hill4(x,xdata);

	xdata = sort(xdata);

	lh = [lh plot(xdata,hill4(x,xdata),'Color',c(i,:))];
	r2 = rsquare(LNFit,ydata);
	L = [L oval(r2,2)];
end
xlabel('Filtered Stimulus (a.u.)')
ylabel('Static Nonlinearity (Hz)')
legend(lh,L,'Location','southeast')

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end



%% Neuron Responses: Fast Adaptation 
% How well do LN models predict the response of neurons in these paradigms? Do we see evidence of fast gain adaptation in these data sets? Now we perform our gain analysis of the prediction using standard methods described elsewhere in this repo. We do the analysis for each stimulus set. 


for i = 1:length(LNModel)
	ph = [];
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	this_orn=mean2(combined_data.fA(:,plot_these));
	this_pid=mean2(combined_data.PID(plot_these,:));
	time = 3e-3*(1:length(this_orn));

	history_lengths = (3*floor(1000*logspace(-2,1,30)/3))/1e3;
	example_history_length = 0.135;

	f2=figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
	ph(3) = subplot(1,2,1); hold on 
	axis square
	ph(4) = subplot(1,2,2); hold on

	% remove trend
	ff = fit(time(:),this_pid(:),'poly2');
	this_pid = this_pid - ff(time)';


	clear x
	x.response = this_orn(a:z-32); % the 32 is to account for the acausal part of the filter
	x.prediction = LNModel(i).LNFit;
	x.stimulus = this_pid(a:z-32);
	x.time = time(a:z-32);
	x.filter_length = 299;

	if redo
		[p_LN,l,h] = GainAnalysis4(x,history_lengths,example_history_length,ph);
		s=abs(l-h);
		s(p_LN(1,:)>0.05)=NaN;
		[~,loc]=max(s);

		% save it for later
		LNModel(i).ehl = history_lengths(loc);
		LNModel(i).p = p_LN;

	else
		GainAnalysis4(x,history_lengths,LNModel(i).ehl,ph,LNModel(i).p);
	end

	xlabel(ph(3),'LN Prediction (Hz)')
	set(ph(4),'XScale','log')
	title(ph(4),paradigm_names{i})

	if being_published

		snapnow;
		delete(f2);
	end
end







%% Version Info
% The file that generated this document is called:
disp(mfilename)

%%
% and its md5 hash is:
Opt.Input = 'file';
disp(DataHash(strcat(mfilename,'.m'),Opt))

%%
% This file should be in this commit:
[status,m]=unix('git rev-parse HEAD');
if ~status
	disp(m)
end
