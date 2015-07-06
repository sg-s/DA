% LFP_MSG.m
% LFP analysis of mean shifted gaussian odour inputs
%
% created by Srinivas Gorur-Shandilya at 10:48 , 03 July 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end
tic

%      ######  ######## #### ##     ## ##     ## ##       ##     ##  ######  
%     ##    ##    ##     ##  ###   ### ##     ## ##       ##     ## ##    ## 
%     ##          ##     ##  #### #### ##     ## ##       ##     ## ##       
%      ######     ##     ##  ## ### ## ##     ## ##       ##     ##  ######  
%           ##    ##     ##  ##     ## ##     ## ##       ##     ##       ## 
%     ##    ##    ##     ##  ##     ## ##     ## ##       ##     ## ##    ## 
%      ######     ##    #### ##     ##  #######  ########  #######   ######  


%% Stimulus Characterisation
% First, we show that we are able to deliver Gaussian-distributed odour stimuli, and that we are able to vary the means of the distributions of these stimuli. 


[PID, LFP, fA, paradigm, orn, AllControlParadigms, paradigm_hashes] = consolidateData('/local-data/DA-paper/LFP-MSG/',1);


% sort the paradigms sensibly
sort_value = [];
for i = 1:length(AllControlParadigms)
	sort_value(i) = mean(mean(AllControlParadigms(i).Outputs));
end

[~,idx] = sort(sort_value);
AllControlParadigms = AllControlParadigms(idx);
paradigm_new = paradigm*NaN;
for i = 1:length(idx)
	paradigm_new(paradigm == i) = idx(i);
end
paradigm = paradigm_new;


% throw out trials where we didn't record the LFP, for whatever reason
not_LFP = find((max(abs(LFP))) < 0.1);
LFP(:,not_LFP) = NaN;

%%
% The following figure shows the distribution of the inputs, for the terminal 40seconds of each 60second presentation. 
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
c = parula(1+length(unique(paradigm)));
for i = 1:length(unique(paradigm))
	hist_this = PID(20e3:end,paradigm==i);
	xx =  linspace(min(min(hist_this)),max(max(hist_this)),50);
	y = NaN(sum(paradigm==i),50);
	for j = 1:sum(paradigm==i)
		y(j,:) = hist(hist_this(:,j),xx);
		y(j,:) = y(j,:)/sum(y(j,:));
	end
	% get everything on the same x axis
	if width(y) > 1
		errorShade(xx,mean(y),std(y)/sqrt(width(y)),'Color',c(i,:),'LineWidth',5);
	else
		plot(xx,y,'Color',c(i,:));
	end
end
xlabel('Stimulus (V)')
ylabel('p(stimulus)')
PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%%
% How reproducible is the stimulus? In the following figure, we show the same stimulus from the entire dataset, over all the trails, and all the flies we have. The shading shows the standard error of the mean.

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot_this = PID(40e3:end,paradigm==1);

time = 1e-3*(1:length(plot_this)) + 40;
errorShade(time,mean2(plot_this),std(plot_this')/sqrt(width(plot_this)));
xlabel('Time (s)')
ylabel('PID (V)')
title(strcat('Stimulus reproducibility: n = ',oval(width(plot_this))))

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

% 								##       ######## ########  
% 								##       ##       ##     ## 
% 								##       ##       ##     ## 
% 								##       ######   ########  
% 								##       ##       ##        
% 								##       ##       ##        
% 								######## ##       ##        


%% Local Field Potential
% We now look at the LFP. Here is an example neuron, showing how the LFP changes with the different paradigms. In the following figure, we plot the raw LFP traces, downsampled to 1kHz from the actual 10kHz trace, and whose baselines (with no odour) have been set to zero. 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
c = parula(1+length(unique(paradigm(:,orn == 7))));
time = 1e-3*(1:length(LFP));
for i = 1:length(c)
	plot_this = find(orn == 7 & paradigm == i);
	for j = 1:length(plot_this)
		this_LFP = LFP(:,plot_this(j));
		this_LFP = this_LFP - mean(this_LFP(1e3:5e3));
		plot(time,this_LFP,'Color',c(i,:))
	end
end

xlabel('Time (s)')
ylabel('LFP (100x V)')
PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%% 
% Some things are clear from this trace: the LFPs look somewhat similar, even for very different odour concentrations. 

%%
% Since we are interested in the response of the neuron to the odour flicker, we only consider the part of the trace corresponding to the odour flicker. Furthermore, we bandpass the LFP trace to remove spikes and to remove these slow fluctuations. This data is from one neuron.

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
these_paradigms = unique(paradigm(:,orn == 7));
c = parula(1+length(these_paradigms));
time = 1e-3*(1:length(LFP(30e3:40e3,1))) + 30;
for i = 1:length(these_paradigms)
	plot_this = find(orn == 7 & paradigm == these_paradigms(i));
	for j = 1:length(plot_this)
		this_LFP = LFP(:,plot_this(j));
		this_LFP = this_LFP(30e3:40e3);
		this_LFP = filter_trace(this_LFP,1000,10);
		plot(time,this_LFP,'Color',c(i,:))
	end
end

xlabel('Time (s)')
ylabel('LFP (100x V)')
PrettyFig;

if being_published
	snapnow
	delete(gcf)
end


%%
% Here we something very strange. THe LFP seems to go 180° out of phase as we move to higher concentrations. What is going on? To be very clear, we plot the stimulus and the LFP for these traces, and normalise them to visualise them togethe: (red is the LFP, and black is the PID).

figure('outerposition',[0 0 1200 700],'PaperUnits','points','PaperSize',[1200 700]); hold on
for i = 1:length(these_paradigms)
	subplot(2,3,i), hold on
	plot_this = find(orn == 7 & paradigm == these_paradigms(i));
	for j = 1:length(plot_this)
		this_LFP = LFP(:,plot_this(j));
		this_LFP = this_LFP(37e3:40e3);
		this_LFP = filter_trace(this_LFP,1000,10);
		this_LFP = this_LFP - mean(this_LFP);
		this_LFP = this_LFP/std(this_LFP);
		time = 1e-3*(1:length(this_LFP)) + 37;
		plot(time,this_LFP,'r')

		this_PID = PID(:,plot_this(j));
		this_PID = this_PID(37e3:40e3);
		this_PID = this_PID - mean(this_PID);
		this_PID = this_PID/std(this_PID);
		plot(time,this_PID,'k')


	end
	title(AllControlParadigms(these_paradigms(i)).Name)
end

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%%
% WTF. 

% ##       ######## ########     ######## #### ##       ######## ######## ########   ######  
% ##       ##       ##     ##    ##        ##  ##          ##    ##       ##     ## ##    ## 
% ##       ##       ##     ##    ##        ##  ##          ##    ##       ##     ## ##       
% ##       ######   ########     ######    ##  ##          ##    ######   ########   ######  
% ##       ##       ##           ##        ##  ##          ##    ##       ##   ##         ## 
% ##       ##       ##           ##        ##  ##          ##    ##       ##    ##  ##    ## 
% ######## ##       ##           ##       #### ########    ##    ######## ##     ##  ######  


%%
% To figure out what's going on, we back out filters from the stimulus to the LFP for each of these cases. In the following figure, we back out filters from all the data we have, and plot them colour coded by paradigm:


K = NaN(1e3,length(orn));
for i = 1:length(orn)
	resp = LFP(30e3:end,i);
	rm_this = isnan(resp);
	resp(rm_this) = [];
	if length(resp) 
		resp = filter_trace(resp,1e3,10);
		stim = PID(30e3:end,i);
		stim(rm_this) = [];
		stim(1:500) = [];
		resp(end-499:end) = [];
		K(:,i) = FitFilter2Data(stim,resp,[],'reg=1;','filter_length=999;');
	end
end


c= parula(max(paradigm)+1);
l = [];
figure('outerposition',[0 0 700 700],'PaperUnits','points','PaperSize',[1200 700]); hold on
for i = 1:max(paradigm)
	time = 1e-3*(1:501)-.2;
	plot_this = find(paradigm == i);
	plot_this = setdiff(plot_this,find(isnan(sum(K))));
	if length(plot_this) > 1
		l(i) = errorShade(time,mean2(K(300:800,plot_this)),std(K(300:800,plot_this)')/length(plot_this),'Color',c(i,:));
	else
		l(i) = plot(time,K(300:800,plot_this),'Color',c(i,:));
	end
end
legend(l,{AllControlParadigms.Name})
xlabel('Lag (s)')
ylabel('PID \rightarrow LFP Filter')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%% 
% A LFP filter with a positive lobe is unheard of. It doesn’t make any sense, and contradicts published data and our own results from other experiments. How widespread is this bizarre positive lobe? In the following figure, we break up the filters for the lowest dose into each neuron, to see if all neurons show the same crazy behaviour: 

c = parula(max(orn)+1);
l = [];
L = {};
figure('outerposition',[0 0 700 700],'PaperUnits','points','PaperSize',[1200 700]); hold on
for i = 1:max(orn)
	time = 1e-3*(1:501)-.2;
	plot_this = find(paradigm == 1 & orn == i);
	plot_this = setdiff(plot_this,find(isnan(sum(K))));
	if length(plot_this) > 1
		l(i) = errorShade(time,mean2(K(300:800,plot_this)),std(K(300:800,plot_this)')/length(plot_this),'Color',c(i,:));
	else
		l(i) = plot(time,K(300:800,plot_this),'Color',c(i,:));
	end
	L{i} = strcat('ORN ',oval(i));
end
legend(l,L)
xlabel('Lag (s)')
ylabel('PID \rightarrow LFP Filter')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%%
% So every ORN shows this, meaning that the problem is either some weird thing with the stimulus, or the flies are fundamentally unsound. 

%% Spiking Filters
% Do the spiking filters also show this weird inversion? To check, we back out spiking filters for all the data and check as before. The following figure shows the filters backed out for some of the lower concentrations:

K_fA = NaN(1e3,length(orn));
for i = 1:length(orn)
	resp = fA(30e3:end,i);
	rm_this = isnan(resp);
	resp(rm_this) = [];
	if length(resp) 
		stim = PID(30e3:end,i);
		stim(rm_this) = [];
		stim(1:500) = [];
		resp(end-499:end) = [];
		K_fA(:,i) = FitFilter2Data(stim,resp,[],'reg=1;','filter_length=999;');
	end
end


c= parula(max(paradigm)+1);
l = [];
figure('outerposition',[0 0 700 700],'PaperUnits','points','PaperSize',[1200 700]); hold on
for i = 1:4
	time = 1e-3*(1:501)-.2;
	plot_this = find(paradigm == i);
	plot_this = setdiff(plot_this,find(isnan(sum(K_fA))));
	if length(plot_this) > 1
		l(i) = errorShade(time,mean2(K_fA(300:800,plot_this)),std(K_fA(300:800,plot_this)')/length(plot_this),'Color',c(i,:));
	elseif length(plot_this) == 1
		l(i) = plot(time,K_fA(300:800,plot_this),'Color',c(i,:));
	else
		l(i) = plot(NaN,NaN);
	end
end
legend(l,{AllControlParadigms.Name})
xlabel('Lag (s)')
ylabel('PID \rightarrow Firing Rate')

PrettyFig;

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
disp(DataHash(strcat(mfilename,'.m'),Opt))

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

