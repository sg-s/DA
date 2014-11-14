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
if redo
	combined_data = ReduceORNData(data_root,allfiles);
end

paradigm_names = unique(combined_data.paradigm);


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


nbins = 50;
histx = [];
histy = [];
paradigm = [];
t_start = 15;
t_stop = 55;
dt = 3e-3;
all_pid = [];

% assemble all histograms
for i = 1:length(paradigm_names)
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	for j = 1:length(plot_these)
		this_pid = combined_data.PID(plot_these(j),:);
		this_pid = this_pid(round(t_start/dt):floor(t_stop/dt));
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
a = floor(15/3e-3);
z = floor(55/3e-3);
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

%% Neuron Responses: Input-Output Curve Changes
% How do ORNs change their input-output curves when presented with these different stimuli? If the ORN is a perfect sensor, that adapts perfectly, it would move and stretch its I/O curve so that it is equal to cumulative density function of the stimulus. However, we know that ORNs don't adapt perfectly, so they must do something else. What is it? If the change dominated by a lateral movement of a stretch?

% this stores the LN model parameters for all the data
clear LNModel
LNModel.K = [];
LNModel.H = []; % stores hill parameters
LNModel.LinearFit = [];
LNModel.LNFit = [];
LNModel.H_domain = [];

a = floor(15/3e-3);
z = floor(55/3e-3);

for i = 1:length(paradigm_names)
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	this_orn=mean2(combined_data.fA(:,plot_these));
	this_pid = mean2(combined_data.PID(plot_these,:));

	time = 3e-3*(1:length(this_orn));

	this_pid = this_pid(a:z);
	this_orn = this_orn(a:z);

	[K,~,filtertime] = FindBestFilter(this_pid,this_orn,[],'filter_length=299;');
	filtertime = filtertime*mean(diff(time));
	K = K/max(K);
	LinearFit = convolve(time,this_pid,K,filtertime) + mean(this_orn);

	xdata = LinearFit(:);
	ydata = this_orn(:);

	ydata(isnan(xdata)) = [];
	xdata(isnan(xdata)) = [];


	fo=optimset('MaxFunEvals',1000,'Display','none');
	x = lsqcurvefit(@hill4,[max(ydata) 2 2 0],xdata,ydata,[max(ydata)/2 2 1 0],[2*max(ydata) max(ydata) 10 max(ydata)],fo);
	LNFit = hill4(x,xdata);

	% save all of this for later
	LNModel(i).K = K;
	LNModel(i).H = x;
	LNModel(i).LinearFit = LinearFit;
	LNModel(i).LNFit = LNFit;
	LNModel(i).H_domain = xdata;
	LNModel(i).H_range =  ydata;
	LNModel(i).LNFit_r2 = rsquare(LNFit,ydata);
	LNModel(i).LinearFit_r2 = rsquare(xdata,ydata);
end


%%
% The following figure shows each of the ORN responses to each stimulus together with the best linear fits:


figure('outerposition',[0 0 1400 600],'PaperUnits','points','PaperSize',[1400 600]); hold on
for i = 1:6
	subplot(2,3,i), hold on
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	this_orn=mean2(combined_data.fA(:,plot_these));
	this_orn = this_orn(a:z);
	time = 3e-3*(1:length(this_orn));
	time = time+a*3e-3;
	plot(time,this_orn,'k');
	time = 3e-3*(1:length(LNModel(i).LinearFit));
	time = time+a*3e-3;
	lh=plot(time,LNModel(i).LinearFit,'r');
	set(gca,'XLim',[30 40],'YLim',[0 55])
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
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	this_orn=mean2(combined_data.fA(:,plot_these));
	this_orn = this_orn(a:z);
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

return
%%
% Are these changes largely subtractive or largely divisive? In the following figure, we attempt to fit the data using only subtractive, only divisive, or mixed models. 

Lh= [];
L = {};
figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
for i = 1:3
	p(i) =  subplot(1,3,i); hold on
	% lowest mean case 
	x=LNModel(i).H_domain;
	y=LNModel(i).H_range;
	scatter(p(i),x,y,'k')
	fo=optimset('MaxFunEvals',1000,'Display','none');
	hf0 = lsqcurvefit(@hill4,[max(y) 2 2 0],x(:),y(:),[max(y)/2 0 1 0],[2*max(y) max(y) 10 10],fo);
	Lh(1) = plot(p(i),sort(x),hill4(hf0,sort(x)),'k');
	L{1} = oval(rsquare(y,hill4(hf,x)),2);
end

backgrounds = [5 6];
title(p(1),'Divisive Only')
for i = 1:length(backgrounds)
	baseline = mean(data(backgrounds(i)).PID(2:70,:));
	x=max(data(backgrounds(i)).PID);
	y=max(data(backgrounds(i)).ORN);
	scatter(p(1),x,y,32,C(i,:))
	fo=optimset('MaxFunEvals',1000,'Display','none');
	hf = lsqcurvefit(@hill4,[max(y) hf0(2) 2 0],x(:),y(:),[max(y)/2 hf0(2) 0 0],[2*max(y) hf0(2)+1e-6 10 10],fo);
	Lh(i+1) = plot(p(1),sort(x),hill4(hf,sort(x)),'Color',C(i,:));
	div_only_hf = hf;
	L{i+1} = oval(rsquare(y,hill4(hf,x)),2);
end
xlabel(p(1),'Stimulus (V)')
ylabel(p(1),'Response (Hz)')
set(p(1),'XScale','log')
legend(p(1),Lh,L,'location','northwest')


title(p(2),'Subtractive Only')
for i = 1:length(backgrounds)
	baseline = mean(data(backgrounds(i)).PID(2:70,:));
	x=max(data(backgrounds(i)).PID);
	y=max(data(backgrounds(i)).ORN);
	scatter(p(2),x,y,32,C(i,:))
	fo=optimset('MaxFunEvals',1000,'Display','none');
	hf = lsqcurvefit(@hill4,[max(y) hf0(2) hf0(3) 0],x(:),y(:),[max(y)/2 0 hf0(3) 0],[2*max(y) 10 hf0(3)+1e-6 10],fo);
	Lh(i+1) = plot(p(2),sort(x),hill4(hf,sort(x)),'Color',C(i,:));
	sub_only_hf = hf;
	L{i+1} = oval(rsquare(y,hill4(hf,x)),2);
end
xlabel(p(2),'Stimulus (V)')
set(p(2),'XScale','log')
legend(p(2),Lh,L,'location','northwest')


title(p(3),'Mixed Model')
for i = 1:length(backgrounds)
	baseline = mean(data(backgrounds(i)).PID(2:70,:));
	x=max(data(backgrounds(i)).PID);
	y=max(data(backgrounds(i)).ORN);
	scatter(p(3),x,y,32,C(i,:))
	fo=optimset('MaxFunEvals',1000,'Display','none');
	hf = lsqcurvefit(@hill4,[max(y) 2 2 0],x(:),y(:),[max(y)/2 0 0 0],[2*max(y) max(y) 10 10],fo);
	Lh(i+1)=plot(p(3),sort(x),hill4(hf,sort(x)),'Color',C(i,:));
	mixed_hf = hf;
	L{i+1} = oval(rsquare(y,hill4(hf,x)),2);
end
xlabel(p(3),'Stimulus (V)')
set(p(3),'XScale','log')
legend(p(3),Lh,L,'location','northwest')



return

%% Neuron Responses: Fast Adaptation 
% How well do LN models predict the response of neurons in these paradigms? Do we see evidence of fast gain adaptation in these data sets? 

%%
% First, we fit a LN model to the data set: 



% build a simple linear model




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
