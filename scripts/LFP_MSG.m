% LFP_MSG.m
% LFP analysis of mean shifted gaussian odour inputs
%
% created by Srinivas Gorur-Shandilya at 10:48 , 03 July 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

% add homebrew path
path1 = getenv('PATH');
path1 = [path1 ':/usr/local/bin'];
setenv('PATH', path1);

calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
		unix(['tag -a publish-failed ',which(mfilename)]);
		unix(['tag -r published ',which(mfilename)]);
	end
end
tic

%% LFP Analysis of Mean Shifted Gaussians
% This document analyses how the LFP varies with Gaussian-distributed stimuli with similar variances but different means. Previously, we showed that the ORN firing rate compresses with increasing mean stimulus, similar to what the Weber-Fechner Law predicts. Here, we do the same for the LFP. 

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

% throw our bad traces
bad_trials = (sum(fA) == 0 | isnan(sum(fA)) |  isnan(sum(LFP)));
LFP(:,bad_trials) = [];
PID(:,bad_trials) = [];
fA(:,bad_trials) = [];
paradigm(bad_trials) = [];
orn(bad_trials) = [];

% band pass all the LFP
filtered_LFP = LFP;
for i = 1:width(LFP)
	filtered_LFP(:,i) = bandPass(LFP(:,i),1000,10);
end


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
% How reproducible is the stimulus? In the following figure, we show the same stimulus from the entire dataset, over all the trails, and all the flies we have. Each trace is stimulus presented to a different ORN. 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot_this = PID(40e3:end,paradigm==1);

time = 1e-3*(1:length(plot_this)) + 40;
plot(time,plot_this,'LineWidth',1,'Color',c(1,:));
xlabel('Time (s)')
ylabel('PID (V)')
title(strcat('Stimulus reproducibility: n = ',oval(width(plot_this))))

prettyFig;

if being_published
	snapnow
	delete(gcf)
end


%    ######## #### ########  #### ##    ##  ######       ######      ###    #### ##    ## 
%    ##        ##  ##     ##  ##  ###   ## ##    ##     ##    ##    ## ##    ##  ###   ## 
%    ##        ##  ##     ##  ##  ####  ## ##           ##         ##   ##   ##  ####  ## 
%    ######    ##  ########   ##  ## ## ## ##   ####    ##   #### ##     ##  ##  ## ## ## 
%    ##        ##  ##   ##    ##  ##  #### ##    ##     ##    ##  #########  ##  ##  #### 
%    ##        ##  ##    ##   ##  ##   ### ##    ##     ##    ##  ##     ##  ##  ##   ### 
%    ##       #### ##     ## #### ##    ##  ######       ######   ##     ## #### ##    ## 


%% Gain in Firing rate
% We now look at how the firing rate gain changes with changes in the mean of the stimulus. We want to check if we can reproduce our earlier results, where we showed that the gain of the firing rates goes down with increasing mean stimulus, and is well fit by a power law with exponent close to -1. 

a = 10e3; z = 50e3;
[K2,fA_pred,fA_gain,fA_gain_err] = extractFilters(PID,fA,'use_cache',true,'a',a,'z',z);


ss = 50;
c = parula(max(paradigm)+1);
l = [];
filtertime = 1e-3*(1:length(K2))-.1;
figure('outerposition',[0 0 1400 500],'PaperUnits','points','PaperSize',[1400 500]); hold on
subplot(1,10,1:5), hold on
for i = 1:max(paradigm)
	plot_this = find(paradigm == i);
	plot_this = setdiff(plot_this,find(isnan(sum(K2))));
	y = mean2(fA(a:z,plot_this));
	x = mean2(fA_pred(a:z,plot_this));
	l(i) = plot(x(1:ss:end),y(1:ss:end),'.','Color',c(i,:));
end

legend(l,{AllControlParadigms.Name},'Location','southeastoutside')
xlabel('Linear Prediction')
ylabel('Firing rate (Hz)')

subplot(1,10,7:10), hold on
for i = 1:width(fA)
	plot(mean(PID(a:z,i)),fA_gain(i),'+','Color',c(paradigm(i),:))
end

x = mean(PID(a:z,:)); x = x(:);
y = fA_gain(:);

xx = NaN(length(unique(paradigm)),1);
yy = NaN(length(unique(paradigm)),1);
ww = NaN(length(unique(paradigm)),1);

for i = 1:length(yy)
	xx(i) = mean(x(paradigm==i));
	yy(i) = mean(y(paradigm==i));
	ww(i) = 1./sem(y(paradigm==i));
end


ff = fit(xx(:),yy(:),'power1','Weights',ww);
clear l
l(1) = plot(sort(xx),ff(sort(xx)),'k--');
L = {};
L{1} = ['y = \alpha (x^\beta) ,\beta = ', oval(ff.b),', r^2=' oval(rsquare(ff(x),y))];


fo = fitoptions('power1');
fo.Upper = [NaN -1];
fo.Lower = [NaN -1];
fo.Weights = ww;
ff = fit(xx(:),yy(:),'power1',fo);
l(2) = plot(sort(xx),ff(sort(xx)),'k');
L{2} = ['y = \alpha (x^\beta) ,\beta := -1, r^2=' oval(rsquare(ff(x),y))];

legend(l,L)
set(gca,'XScale','log','YScale','log')
xlabel('Mean Stimulus (V)')
ylabel('ORN Gain (Hz/V)')

prettyFig;

if being_published
	snapnow
	delete(gcf)
end


%%
% So this is very nice, as we can get the same result all over again, and the exponent of the fit is also very close to the old value, and very close to the theoretical prediction (-1).  


% 								##       ######## ########  
% 								##       ##       ##     ## 
% 								##       ##       ##     ## 
% 								##       ######   ########  
% 								##       ##       ##        
% 								##       ##       ##        
% 								######## ##       ##        


%% Local Field Potential
% We now look at the LFP. Here is an example neuron, showing how the LFP changes with the different paradigms. In the following figure, we plot the raw LFP traces, downsampled to 1kHz from the actual 10kHz trace, and whose baselines (with no odour) have been set to zero. This figure shows all the features of the phenomenology we are interested in: 
% 
% # The maximum absolute value of the LFP scales with the magnitude of the stimulus applied (yellow traces are lower than blue traces)
% # LFP gain at high stimuli is lower than LFP gain at low stimuli (standard deviation of yellow traces is less than that of blue traces)
% # Gain controls happens over the course of the experiment. At t=10s, the yellow trace still responds to the Gaussian distributed flicker, but at t>40s, those fluctuations have been damped out. 

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
		plot(time,filtered_LFP(10e3:55e3,plot_this(j)),'Color',c(i,:))
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
% To figure out what's going on, we back out filters from the stimulus to the LFP for each of these cases, make linear predictions, and then estimate gain, just like we did with the firing rates. 

a = 10e3; z = 50e3;
[K1,LFP_pred,LFP_gain,LFP_gain_err] = extractFilters(PID,filtered_LFP,'use_cache',true,'a',a,'z',z);

% how good are these fits?
r2 = NaN(width(LFP),1);
for i = 1:width(LFP)
	r2(i) = rsquare(LFP_pred(a:z,i),filtered_LFP(a:z,i));
end

ss  =50;
c= parula(max(paradigm)+1);
l = [];
filtertime = 1e-3*(1:length(K1))-.1;
figure('outerposition',[0 0 1400 500],'PaperUnits','points','PaperSize',[1400 500]); hold on
subplot(1,10,1:5), hold on
for i = 1:max(paradigm)
	plot_this = find(paradigm == i);
	plot_this = setdiff(plot_this,find(isnan(sum(K1))));
	y = filtered_LFP(a:z,plot_this);
	y = mean2(y);
	x = mean2(LFP_pred(a:z,plot_this));
	
	l(i) = plot(x(1:ss:end),y(1:ss:end),'.','Color',c(i,:));
end
legend(l,{AllControlParadigms.Name},'Location','southeastoutside')
xlabel('Linear Prediction')
ylabel('\DeltaLFP (mV)')

subplot(1,10,7:10), hold on
for i = 1:width(LFP)
	plot(mean(PID(a:z,i)),LFP_gain(i),'+','Color',c(paradigm(i),:))
end

x = mean(PID(a:z,:)); x = x(:);
y = LFP_gain(:);

xx = NaN(length(unique(paradigm)),1);
yy = NaN(length(unique(paradigm)),1);
ww = NaN(length(unique(paradigm)),1);

for i = 1:length(yy)
	xx(i) = mean(x(paradigm==i));
	yy(i) = mean(y(paradigm==i));
	ww(i) = 1./sem(y(paradigm==i));
end


ff = fit(xx(:),yy(:),'power1','Weights',ww);
clear l
l(1) = plot(sort(xx),ff(sort(xx)),'k--');
L = {};
L{1} = ['y = \alpha (x^\beta) ,\beta = ', oval(ff.b),', r^2=' oval(rsquare(ff(x),y))];


fo = fitoptions('power1');
fo.Upper = [NaN -1];
fo.Lower = [NaN -1];
fo.Weights = ww;
ff = fit(xx(:),yy(:),'power1',fo);
l(2) = plot(sort(xx),ff(sort(xx)),'k');
L{2} = ['y = \alpha (x^\beta) ,\beta := -1, r^2=' oval(rsquare(ff(x),y))];

legend(l,L)
set(gca,'XScale','log','YScale','log')
xlabel('Mean Stimulus (V)')
ylabel('LFP Gain')

prettyFig;

if being_published
	snapnow
	delete(gcf)
end


%% LFP to Firing Analysis
% In this section we analyse the transformation from the LFP to the firing, and see how gain here varies with the stimulus. 

a = 10e3; z = 50e3;
[K3,fA_pred,fA_gain,fA_gain_err] = extractFilters(-LFP,fA,'band_pass_x',true,'use_cache',true,'a',a,'z',z);

% how good are these fits?
r2 = NaN(width(LFP),1);
for i = 1:width(LFP)
	r2(i) = rsquare(fA_pred(a:z,i),fA(a:z,i));
end

ss  =50;
c= parula(max(paradigm)+1);
l = [];
filtertime = 1e-3*(1:length(K3))-.1;

figure('outerposition',[0 0 1400 500],'PaperUnits','points','PaperSize',[1400 500]); hold on
subplot(1,10,1:5), hold on
for i = 1:max(paradigm)
	plot_this = find(paradigm == i);
	plot_this = setdiff(plot_this,find(isnan(sum(K3))));
	y = fA(a:z,plot_this);
	y = mean2(y);
	x = mean2(fA_pred(a:z,plot_this));
	l(i) = plot(x(1:ss:end),y(1:ss:end),'.','Color',c(i,:));
end
legend(l,{AllControlParadigms.Name},'Location','southeastoutside')
xlabel('-Linear Prediction')
ylabel('Firing Rate (Hz)')

subplot(1,10,7:10), hold on
for i = 1:width(LFP)
	plot(mean(PID(a:z,i)),abs(fA_gain(i)),'+','Color',c(paradigm(i),:))
end

x = mean(PID(a:z,:)); x = x(:);
y = fA_gain(:);

xx = NaN(length(unique(paradigm)),1);
yy = NaN(length(unique(paradigm)),1);
ww = NaN(length(unique(paradigm)),1);

for i = 1:length(yy)
	xx(i) = mean(x(paradigm==i));
	yy(i) = mean(y(paradigm==i));
	ww(i) = 1./sem(y(paradigm==i));
end


ff = fit(xx(:),yy(:),'power1','Weights',ww);
clear l
l(1) = plot(sort(xx),ff(sort(xx)),'k--');
L = {};
L{1} = ['y = \alpha (x^\beta) ,\beta = ', oval(ff.b) ,' , r^2=' oval(rsquare(ff(x),y))];

legend(l,L,'Location','northwest')
set(gca,'XScale','log','YScale','log')
xlabel('Mean Stimulus (V)')
ylabel('LFP -> Firing Gain')
set(gca,'YLim',[50 1000])

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
disp(dataHash(strcat(mfilename,'.m'),Opt))

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

% tag the file as being published 
% add homebrew path
path1 = getenv('PATH');
path1 = [path1 ':/usr/local/bin'];
setenv('PATH', path1);

if being_published
	unix(['tag -a published ',which(mfilename)]);
	unix(['tag -r publish-failed ',which(mfilename)]);
end

