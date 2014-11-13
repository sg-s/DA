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

%% Mean Shifted Gaussians
% How to ORNs respond to mean shifted gaussians? Specifcally, how do they vary their input-output curve? Is the adaptation to this mean optimal (a la Laughlin etc)? Can we find evidence for fast gain adaptation in the curves themselves? 

data_root = '/local-data/DA-paper/fast-flicker/orn/';
allfiles  = dir(strcat(data_root,'*.mat'));

%% Stimulus Distributions 
% Three gaussians with different means are chosen. The following figure shows the distributions for every trial of the data analysed here, colour coded by the experimental paradigm. 


nbins = 50;
histx = [];
histy = [];
paradigm = [];
t_start = 15;
t_stop = 55;
dt = 1e-4;
all_pid = [];

% assemble all histograms
for i = 1:length(allfiles)
	load(strcat(data_root,allfiles(i).name))
	for j = 1:length(data)
		if ~isempty(data(j).PID)
			if length(data(j).PID) > t_stop/dt
				for k = 1:width(data(j).PID)
					this_pid = data(j).PID(k,floor(t_start/dt):floor(t_stop/dt));

					[y,x] = hist(this_pid,nbins);
					if mean(x) > .5
						histx = [histx; x];
						histy = [histy; y];
						paradigm = [paradigm j];
						all_pid = [all_pid; this_pid];
						
					end
					
				end
			end
		end
	end
end



c = [0 0 1; 0 1 0; 1 0 0];
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

c = [0 0 1; 0 1 0; 1 0 0];
paradigm = paradigm - min(paradigm);
paradigm = paradigm + 1;
time = 1:length(all_pid);
time = time*dt;

figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
for i = 1:length(paradigm)
	plot(time(1:50:end),all_pid(i,1:50:end),'Color',c(paradigm(i),:))
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

% combine all data
if redo
	combined_data = ReduceORNData(data_root,allfiles);
end


%% Neuron Responses: Overview
% The following figure shows the responses of the ORNs to this stimuli, and their distributions. 

figure('outerposition',[0 0 1400 500],'PaperUnits','points','PaperSize',[1400 500]); hold on
subplot(1,5,1:4), hold on
for i = unique(combined_data.paradigm)
	this_orn=mean2(combined_data.fA(:,combined_data.paradigm==i));
	time = 3e-3*(1:length(this_orn));
	plot(time,this_orn)
end
ylim = get(gca,'YLim');
a = floor(15/3e-3);
z = floor(55/3e-3);
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')

subplot(1,5,5), hold on
for i = unique(combined_data.paradigm)
	this_orn=mean2(combined_data.fA(:,combined_data.paradigm==i));
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


%% Neuron Responses: Fast Adaptation 
% How well do LN models predict the response of neurons in these paradigms? Do we see evidence of fast gain adaptation in these data sets? 

%%
% First, we fit a LN model to the data set: 

this_pid=mean2(combined_data.PID(combined_data.paradigm==2,:));
this_orn=mean2(combined_data.fA(:,combined_data.paradigm==2));
time = 3e-3*(1:length(this_orn));
a = floor(15/3e-3);
z = floor(55/3e-3);

% build a simple linear model
[K,~,filtertime] = FindBestFilter(this_pid(a:z),this_orn(a:z),[],'filter_length=299;');
filtertime = filtertime*mean(diff(time));
K = K/max(K);
LinearFit = convolve(time,this_pid,K,filtertime);
LinearFit = LinearFit + mean(this_orn(a:z));

xdata = LinearFit(a:z);
ydata = this_orn(a:z);

% crop it to lose NaNs
ydata(isnan(xdata)) = [];
xdata(isnan(xdata)) = [];

xdata = xdata(:);
ydata = ydata(:);

fo=optimset('MaxFunEvals',1000,'Display','none');
x = lsqcurvefit(@hill,[max(ydata) 2 2],xdata,ydata,[max(ydata)/2 2 1],[2*max(ydata) max(ydata) 10],fo);
% save this for later
LN_pred = hill(x,LinearFit);


figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
subplot(2,2,1), hold on
plot(filtertime,K,'k')
xlabel('Filter Lag (s)')
ylabel('Filter Amplitude (norm)')

subplot(2,2,2), hold on
plot(xdata,ydata,'.','Color',[.9 .9 .9])
plot(sort(xdata),hill(x,sort(xdata)),'k')
xlabel('Filter Output (a.u.)')
ylabel('ORN Response (Hz)')

subplot(2,2,3:4), hold on
plot(time(a:z),this_orn(a:z),'k')
plot(time(a:z),LN_pred(a:z),'r')
set(gca,'XLim',[30 40])

xlabel('Time (s)')
ylabel('Firing rate (Hz)')
legend({'Data','LN Prediction'})

PrettyFig;
if being_published
	snapnow;
	delete(gcf);
end

%% 
% Now we perform our gain analysis of the prediction using standard methods described elsewhere in this repo.
ph = [];
this_pid=mean2(combined_data.PID(combined_data.paradigm==2,:));
this_orn=mean2(combined_data.fA(:,combined_data.paradigm==2));
time = 3e-3*(1:length(this_orn));

history_lengths = (3*floor(1000*logspace(-2,1,30)/3))/1e3;
example_history_length = 0.135;

f2=figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ph(3) = subplot(1,2,1); hold on 
axis square
ph(4) = subplot(1,2,2); hold on


s = floor(15/3e-3);
z = floor(55/3e-3);

clear x
x.response = this_orn(s:z);
x.prediction = LN_pred(s:z);
x.stimulus = this_pid(s:z);
x.time = time(s:z);
x.filter_length = 299;

if redo
	[p_LN,l,h] = GainAnalysis4(x,history_lengths,example_history_length,ph);
	s=abs(l-h);
	s(p_LN(1,:)>0.05)=NaN;
	[~,loc]=max(s);
	example_history_length_LN = history_lengths(loc);

else
	GainAnalysis4(x,history_lengths,example_history_length_LN,ph,p_LN);
end

xlabel(ph(3),'LN Prediction (Hz)')
set(ph(4),'XScale','log')

if being_published
	snapnow;
	delete(f1);

	snapnow;
	delete(f2);
end



%% Neuron Responses: Step Adaptation 
% This dataset also allows us to look at how ORNs respond to a step on of stimulus in great detail. 


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
