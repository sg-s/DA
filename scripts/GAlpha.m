% GAlpha.m
% 
% created by Srinivas Gorur-Shandilya at 8:04 , 26 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


% add homebrew path
path1 = getenv('PATH');
if isempty(strfind(path1,':/usr/local/bin'))
    path1 = [path1 ':/usr/local/bin'];
end
setenv('PATH', path1);

% this code determines if this function is being called
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


%% Comparison with old experiments
% In this section we show the stimulus and the LFP and compare them to an older experiment where we did the same thing with CS flies. The following figure shows that the stimulus is approximately the same, but the LFP are quite different:

[alldata(1).PID, alldata(1).LFP, alldata(1).fA, alldata(1).paradigm, alldata(1).orn, alldata(1).fly, alldata(1).AllControlParadigms, alldata(1).paradigm_hashes] = consolidateData('/local-data/DA-paper/g-alpha/rnai',1);
alldata(1).genotype = 'UAS-g-alpha-RNAi-22a-gal4';

[alldata(2).PID, alldata(2).LFP, alldata(2).fA, alldata(2).paradigm, alldata(2).orn, alldata(2).fly, alldata(2).AllControlParadigms, alldata(2).paradigm_hashes] = consolidateData('/local-data/DA-paper/LFP-MSG/september',1);
alldata(2).genotype = 'CS';


for ai = 1:length(alldata)
	% remove baseline from all PIDs
	for i = 1:width(alldata(ai).PID)
		alldata(ai).PID(:,i) = alldata(ai).PID(:,i) - mean(alldata(ai).PID(1:5e3,i));
	end

	% remove baseline from all LFPs
	for i = 1:width(alldata(ai).LFP)
		alldata(ai).LFP(:,i) = alldata(ai).LFP(:,i) - mean(alldata(ai).LFP(1:5e3,i));
	end

	% sort the paradigms sensibly
	sort_value = [];
	for i = 1:length(alldata(ai).AllControlParadigms)
		sort_value(i) = (mean(alldata(ai).AllControlParadigms(i).Outputs(1,:)));
	end

	[~,idx] = sort(sort_value);


	alldata(ai).AllControlParadigms = alldata(ai).AllControlParadigms(idx);
	alldata(ai).paradigm_hashes = alldata(ai).paradigm_hashes(idx);
	paradigm_new = alldata(ai).paradigm*NaN;
	for i = 1:length(idx)
		paradigm_new(alldata(ai).paradigm == idx(i)) = i;
	end
	alldata(ai).paradigm = paradigm_new;

	% remove "Flicker" from paradigm names
	for i = 1:length(alldata(ai).AllControlParadigms)
		alldata(ai).AllControlParadigms(i).Name = strrep(alldata(ai).AllControlParadigms(i).Name,'Flicker-','');
	end


	% throw out trials where we didn't record the LFP, for whatever reason
	not_LFP = find((max(abs(alldata(ai).LFP))) < 0.1);
	alldata(ai).LFP(:,not_LFP) = NaN;

	% throw out trials where i think we lost the neuron
	lost_neuron = ones(width(alldata(ai).LFP),1);
	for i = 1:width(alldata(ai).LFP)
		temp = alldata(ai).LFP(50e3:60e3,i);
		lost_neuron(i) = nanstd(temp);
	end
	alldata(ai).LFP(:,lost_neuron<.1) = NaN;


	% throw our bad traces
	bad_trials = isnan(sum(alldata(ai).LFP));
	alldata(ai).LFP(:,bad_trials) = [];
	alldata(ai).PID(:,bad_trials) = [];
	alldata(ai).fA(:,bad_trials) = [];
	alldata(ai).paradigm(bad_trials) = [];
	alldata(ai).orn(bad_trials) = [];

	% band pass all the LFP
	alldata(ai).filtered_LFP = alldata(ai).LFP;
	for i = 1:width(alldata(ai).LFP)
		alldata(ai).filtered_LFP(:,i) = bandPass(alldata(ai).LFP(:,i),1000,10);
	end
end

all_paradigm_hashes = unique([alldata.paradigm_hashes]');
c = parula(1+length(all_paradigm_hashes));

figure('outerposition',[0 0 1200 900],'PaperUnits','points','PaperSize',[1200 900]); hold on
subplot(2,2,1), hold on
ai = 1; ss = 10;
for i = 1:length(all_paradigm_hashes)
	this_hash = all_paradigm_hashes{i};
	this_paradigm = find(strcmp(alldata(ai).paradigm_hashes,this_hash));
	plot_this = alldata(ai).PID(1:ss:end,alldata(ai).paradigm==this_paradigm);
	time = 1e-2*(1:length(plot_this));
	plot(time,mean(plot_this,2),'Color',c(this_paradigm,:),'LineWidth',2);
end
ylabel('Stimulus (V)')
xlabel('Time (s)')
set(gca,'YLim',[0 2.5])
title('UAS-G_{\alpha-S}-RNAi;Or22a-GAL4')

subplot(2,2,2), hold on
ai = 2;
for i = 1:length(all_paradigm_hashes)
	this_hash = all_paradigm_hashes{i};
	this_paradigm = find(strcmp(alldata(ai).paradigm_hashes,this_hash));
	if ~isempty(this_paradigm)
		plot_this = alldata(ai).PID(1:ss:end,alldata(ai).paradigm==this_paradigm);
		time = 1e-2*(1:length(plot_this));
		plot(time,mean(plot_this,2),'Color',c(this_paradigm,:),'LineWidth',2);
	end
end
ylabel('Stimulus (V)')
xlabel('Time (s)')
set(gca,'YLim',[0 2.5])
title('Canton S')

subplot(2,2,3), hold on
ai = 1;
for i = 1:length(all_paradigm_hashes)
	this_hash = all_paradigm_hashes{i};
	this_paradigm = find(strcmp(alldata(ai).paradigm_hashes,this_hash));
	if ~isempty(this_paradigm)
		plot_this = alldata(ai).LFP(1:ss:end,alldata(ai).paradigm==this_paradigm);
		time = 1e-2*(1:length(plot_this));
		plot(time,mean(plot_this,2),'Color',c(this_paradigm,:),'LineWidth',2);
	end
end
ylabel('LFP (mV)')
xlabel('Time (s)')
set(gca,'YLim',[-4 1])

subplot(2,2,4), hold on
ai = 2;
for i = 1:length(all_paradigm_hashes)
	this_hash = all_paradigm_hashes{i};
	this_paradigm = find(strcmp(alldata(ai).paradigm_hashes,this_hash));
	if ~isempty(this_paradigm)
		plot_this = alldata(ai).LFP(1:ss:end,alldata(ai).paradigm==this_paradigm);
		time = 1e-2*(1:length(plot_this));
		plot(time,mean(plot_this,2),'Color',c(this_paradigm,:),'LineWidth',2);
	end
end
ylabel('LFP (mV)')
xlabel('Time (s)')
set(gca,'YLim',[-4 1])

prettyFig('fs=14;');

if being_published
	snapnow
	delete(gcf)
end

%%
% It looks like the RNAi-knockdown is more sensitive (so the LFP deflections are bigger), but also that a slow timescale in the CS (the rate of adaptation to the overall step) becomes much slower (i.e., see that the lowest dose elicits a almost perfect adaptation in the CS, but there is hardly any step adaptation in the RNAi mutant). 

%%
% To see this more clearly, we compare the LFP responses paradigm by paradigm. 

figure('outerposition',[0 0 1400 900],'PaperUnits','points','PaperSize',[1400 900]); hold on
common_paradigms = intersect(alldata(1).paradigm_hashes,alldata(2).paradigm_hashes);
% sort by increasing order
n = length(common_paradigms);
conc = zeros(n,1);
for i = 1:n
	this_hash = common_paradigms{i};
	this_paradigm = find(strcmp(alldata(1).paradigm_hashes,this_hash));
	conc(i) = this_paradigm;
end
[~,idx] = sort(conc);
common_paradigms = common_paradigms(idx);

for i = 1:n
	autoPlot(n,i); hold on
	this_hash = common_paradigms{i};
	this_paradigm = find(strcmp(alldata(2).paradigm_hashes,this_hash));
	plot_this = alldata(2).LFP(1:ss:end,alldata(2).paradigm==this_paradigm);
	time = 1e-2*(1:length(plot_this));
	plot(time,mean(plot_this,2),'Color','k','LineWidth',2);

	this_paradigm = find(strcmp(alldata(1).paradigm_hashes,this_hash));
	plot_this = alldata(1).LFP(1:ss:end,alldata(1).paradigm==this_paradigm);
	time = 1e-2*(1:length(plot_this));
	plot(time,mean(plot_this,2),'Color','r','LineWidth',2);

	title(alldata(1).AllControlParadigms(this_paradigm).Name)
end

prettyFig('fs=14;');

if being_published
	snapnow
	delete(gcf)
end

%% LFP Gain
% In this section we analyse the gain in the LFP and compare it to the wild type. 


a = 10e3; z = 50e3;
for ai = 1:length(alldata)
	[alldata(ai).K1,alldata(ai).LFP_pred,alldata(ai).LFP_gain,alldata(ai).LFP_gain_err] = extractFilters(alldata(ai).PID,alldata(ai).filtered_LFP,'use_cache',true,'a',a,'z',z);
end


figure('outerposition',[0 0 600 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
ai = 1;
clear l
for i = 1:width(alldata(ai).LFP)
	l(ai) = plot(mean(alldata(ai).PID(a:z,i)),alldata(ai).LFP_gain(i),'+','Color','r');
end

ai = 2;
for i = 1:width(alldata(ai).LFP)
	l(ai) = plot(mean(alldata(ai).PID(a:z,i)),alldata(ai).LFP_gain(i),'+','Color','k');
end

set(gca,'XScale','log','YScale','log')
xlabel('Mean Stimulus (V)')
ylabel('LFP Gain (mV/V)')
legend(l,{'UAS-g-alpha-RNAi-22a-gal4','CS'})
set(gca,'YLim',[.1 5])

prettyFig;

if being_published
	snapnow
	delete(gcf)
end

%% Kinetics 
% In this section, we check if there are any kinetics changes in the mutant, and compare them to the WT. 

a = 30e3; z = 50e3;
for ai = 1:length(alldata)
	alldata(ai).tau = NaN(width(alldata(ai).PID),1);
	for i = 1:width(alldata(ai).PID)
		x = alldata(ai).PID(a:z,i);
		x = x - mean(x); x = x/std(x);
		y = alldata(ai).filtered_LFP(a:z,i);
		y = y - mean(y); y = y/std(y);
		y = -y;
		temp = xcorr(x,y);
		[~,loc] = max(temp);
		alldata(ai).tau(i) = 20e3 - loc;
	end
end

figure('outerposition',[0 0 600 500],'PaperUnits','points','PaperSize',[600 500]); hold on
ai = 1;
for i = 1:max(alldata(ai).paradigm)
	mean_stim = mean(mean(alldata(ai).PID(a:z,alldata(ai).paradigm==i)));
	mean_tau = mean(alldata(ai).tau(alldata(ai).paradigm==i));
	tau_err = sem(alldata(ai).tau(alldata(ai).paradigm==i));
	errorbar(mean_stim,mean_tau,tau_err,'r')
end

ai = 2;
for i = 1:max(alldata(ai).paradigm)
	mean_stim = mean(mean(alldata(ai).PID(a:z,alldata(ai).paradigm==i)));
	mean_tau = mean(alldata(ai).tau(alldata(ai).paradigm==i));
	tau_err = sem(alldata(ai).tau(alldata(ai).paradigm==i));
	errorbar(mean_stim,mean_tau,tau_err,'k')
end

xlabel('Mean Stimulus (V)')
ylabel('Response delay (ms)')

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

if being_published
	unix(['tag -a published ',which(mfilename)]);
	unix(['tag -r publish-failed ',which(mfilename)]);
end

