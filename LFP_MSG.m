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
% First, we show that we are able to deliver Gaussian-distributed odour stimuli, and that we are able to vary the means of the distributions of these stimuli. The following figure shows the mean of each stimulus drawn from each distribution (left), and the actual distributions themselves (right).  

[PID, LFP, fA, paradigm, orn, fly, AllControlParadigms, paradigm_hashes] = consolidateData('/local-data/DA-paper/LFP-MSG/september',1);


% remove baseline from all PIDs
for i = 1:width(PID)
	PID(:,i) = PID(:,i) - mean(PID(1:5e3,i));
end

% sort the paradigms sensibly
sort_value = [];
for i = 1:length(AllControlParadigms)
	sort_value(i) = (mean(AllControlParadigms(i).Outputs(1,:)));
end

[~,idx] = sort(sort_value);


AllControlParadigms = AllControlParadigms(idx);
paradigm_new = paradigm*NaN;
for i = 1:length(idx)
	paradigm_new(paradigm == idx(i)) = i;
end
paradigm = paradigm_new;

% remove "Flicker" from paradigm names
for i = 1:length(AllControlParadigms)
	AllControlParadigms(i).Name = strrep(AllControlParadigms(i).Name,'Flicker-','');
end


% throw out trials where we didn't record the LFP, for whatever reason
not_LFP = find((max(abs(LFP))) < 0.1);
LFP(:,not_LFP) = NaN;

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,4,1:3), hold on
c = parula(1+length(unique(paradigm)));
for i = 1:length(unique(paradigm))
	plot_this = PID(40e3:45e3,paradigm==i);
	time = 40+1e-3*(1:length(plot_this));
	errorShade(time,mean(plot_this,2),sem(plot_this),'Color',c(i,:),'LineWidth',2);
end
ylabel('Stimulus (V)')
xlabel('Time (s)')
set(gca,'YLim',[0 2])

subplot(1,4,4), hold on
c = parula(1+length(unique(paradigm)));
for i = 1:length(unique(paradigm))
	hist_this = PID(20e3:55e3,paradigm==i);
	xx =  linspace(min(min(hist_this)),max(max(hist_this)),50);
	y = NaN(sum(paradigm==i),50);
	for j = 1:sum(paradigm==i)
		y(j,:) = hist(hist_this(:,j),xx);
		y(j,:) = y(j,:)/sum(y(j,:));
	end
	plot(mean2(y),(xx),'Color',c(i,:));
end

xlabel('p(stimulus)')
set(gca,'YLim',[0 2])
prettyFig;

if being_published
	snapnow
	delete(gcf)
end


%%
% How reproducible is the stimulus? In the following figure, we show the same stimulus from the entire dataset, over all the trails, and all the flies we have. The shading shows the standard error of the mean.

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot_this = PID(40e3:end,paradigm==1);

time = 1e-3*(1:length(plot_this)) + 40;
errorShade(time,mean2(plot_this),std(plot_this')/sqrt(width(plot_this)),'LineWidth',1,'Color',c(1,:));
xlabel('Time (s)')
ylabel('PID (V)')
title(strcat('Stimulus reproducibility: n = ',oval(width(plot_this))))

prettyFig;

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

example_orn = 4;
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
c = parula(1+length(unique(paradigm(:,orn == example_orn))));
time = 1e-3*(1:length(LFP));
for i = 1:length(c)
	plot_this = find(orn == example_orn & paradigm == i);
	for j = 1:length(plot_this)
		this_LFP = LFP(:,plot_this(j));
		this_LFP = this_LFP - mean(this_LFP(1e3:5e3));
		plot(time,this_LFP,'Color',c(i,:))
	end
end

xlabel('Time (s)')
ylabel('LFP (100x V)')
prettyFig;

if being_published
	snapnow
	delete(gcf)
end


%%
% Since we are interested in the response of the neuron to the odour flicker, we only consider the part of the trace corresponding to the odour flicker. Furthermore, we bandpass the LFP trace to remove spikes and to remove these slow fluctuations. This data is from one neuron.

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
these_paradigms = unique(paradigm(:,orn == example_orn));
c = parula(1+length(these_paradigms));
time = 1e-3*(1:length(LFP(10e3:55e3,1))) + 10;
for i = 1:length(these_paradigms)
	plot_this = find(orn == example_orn & paradigm == these_paradigms(i));
	for j = 1:length(plot_this)
		this_LFP = LFP(:,plot_this(j));
		this_LFP = this_LFP(10e3:55e3);
		this_LFP = bandPass(this_LFP,1000,10);
		plot(time,this_LFP,'Color',c(i,:))
	end
end

xlabel('Time (s)')
ylabel('LFP (100x V)')
prettyFig;

if being_published
	snapnow
	delete(gcf)
end



% ##       ######## ########     ######## #### ##       ######## ######## ########   ######  
% ##       ##       ##     ##    ##        ##  ##          ##    ##       ##     ## ##    ## 
% ##       ##       ##     ##    ##        ##  ##          ##    ##       ##     ## ##       
% ##       ######   ########     ######    ##  ##          ##    ######   ########   ######  
% ##       ##       ##           ##        ##  ##          ##    ##       ##   ##         ## 
% ##       ##       ##           ##        ##  ##          ##    ##       ##    ##  ##    ## 
% ######## ##       ##           ##       #### ########    ##    ######## ##     ##  ######  


%%
% To figure out what's going on, we back out filters from the stimulus to the LFP for each of these cases. In the following figure, we back out filters from all the data we have, and plot them colour coded by paradigm:

K = cache(dataHash([PID LFP]));
if isempty(K)
	K = NaN(700,length(orn));
	for i = 1:length(orn)
		textbar(i,length(orn))
		resp = LFP(10e3:55e3,i);
		rm_this = isnan(resp);
		resp(rm_this) = [];
		if length(resp) 
			try
				resp = bandPass(resp,1e3,10);
				stim = PID(10e3:55e3,i);
				stim(rm_this) = [];
				stim(1:300) = [];
				resp(end-299:end) = [];
				[temp,filtertime] = fitFilter2Data(stim,resp,'reg',1,'filter_length',999);
				K(:,i) = temp(201:900);
				filtertime = filtertime(201:900);
			catch

			end
		end
	end
	cache(dataHash(LFP),K);
end

c= parula(max(paradigm)+1);
l = [];
filtertime = 1e-3*(1:length(K))-.1;
figure('outerposition',[0 0 1400 500],'PaperUnits','points','PaperSize',[1400 500]); hold on
subplot(1,10,1:5), hold on
for i = 1:max(paradigm)

	plot_this = find(paradigm == i);
	plot_this = setdiff(plot_this,find(isnan(sum(K))));
	if length(plot_this) > 1
		l(i) = errorShade(filtertime,mean2(K(:,plot_this)),std(K(:,plot_this)')/length(plot_this),'Color',c(i,:));
	else
		l(i) = plot(time,K(:,plot_this),'Color',c(i,:));
	end
end
legend(l,{AllControlParadigms.Name},'Location','southeastoutside')
xlabel('Lag (s)')
ylabel('PID \rightarrow LFP Filter')

subplot(1,10,7:10), hold on

% make linear predictions 
LFP_pred = NaN*LFP;
LFP_gain = NaN*orn;
LFP_gain_err = NaN*orn;
a = 10e3;
z = 55e3;
time = 1e-3*(1:length(LFP));
for i = 1:width(LFP)
	LFP_pred(:,i) = convolve(time,PID(:,i),K(:,i),filtertime);
	% fit lines to estimate gains
	x = LFP_pred(a:z,i);
	y = bandPass(LFP(a:z,i),1e3,10);
	try
		[ff,gof] = fit(x(:),y(:),'poly1');
		LFP_gain_err(i) = 1 - gof.rsquare;
		LFP_gain(i) = ff.p1;
	catch
	end
end

for i = 1:width(LFP)
	plot(mean(PID(a:z,i)),LFP_gain(i),'+','Color',c(paradigm(i),:))
end

% fit a power law to this
x = mean(PID(a:z,:));
y = LFP_gain;
rm_this = isnan(x) | isnan(LFP_gain);
x(rm_this) = []; y(rm_this) = [];
fo = fitoptions('rat01');
fo.StartPoint = [.08 -.17];
ff = fit(x(:),y(:),'power1');
l = plot(sort(x),ff(sort(x)),'k--');
legend(l,['y = \alpha (x^\beta) , r^2=' oval(rsquare(ff(x),y))])

set(gca,'XScale','log','YScale','log')
xlabel('Mean Stimulus (V)')
ylabel('LFP Gain')

prettyFig;

if being_published
	snapnow
	delete(gcf)
end

return

%% Spiking Filters
% We now extract filters for the neuron spiking. 

K2 = cache(dataHash(fA));
if isempty(K2)
	K2 = NaN(700,length(orn));
	for i = 1:length(orn)
		textbar(i,length(orn))
		resp = fA(10e3:55e3,i);
		rm_this = isnan(resp);
		resp(rm_this) = [];
		if length(resp) 
			try
				stim = PID(10e3:55e3,i);
				stim(rm_this) = [];
				stim(1:300) = [];
				resp(end-299:end) = [];
				[temp,filtertime] = fitFilter2Data(stim,resp,'reg',1,'filter_length',999);
				K2(:,i) = temp(201:900);
				filtertime = filtertime(201:900);
			catch

			end
		end
	end
	cache(dataHash(fA),K2);
end


c= parula(max(paradigm)+1);
l = [];
filtertime = 1e-3*(1:length(K2))-.1;
figure('outerposition',[0 0 700 700],'PaperUnits','points','PaperSize',[1200 700]); hold on
for i = 1:max(paradigm)

	plot_this = find(paradigm == i);
	plot_this = setdiff(plot_this,find(isnan(sum(K2))));
	if length(plot_this) > 1
		l(i) = errorShade(filtertime,mean2(K2(:,plot_this)),std(K2(:,plot_this)')/length(plot_this),'Color',c(i,:));
	else
		l(i) = plot(time,K2(:,plot_this),'Color',c(i,:));
	end
end
legend(l,{AllControlParadigms.Name},'Location','southeast')
xlabel('Lag (s)')
ylabel('PID \rightarrow Firing Filter')

prettyFig;

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

