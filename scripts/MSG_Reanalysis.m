% MSG_Reanalysis.m
% 
% created by Srinivas Gorur-Shandilya at 7:41 , 10 November 2015. Contact me at http://srinivas.gs/contact/
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

%% Reanalysis of Mean Shifted Gaussians Experiment 
% In this document we completely re-analyse the mean-shifted Gaussians experiment. The re-analysis is total: I sorted the spikes using new t-SNE based methods (in the old analysis, Mahmut sorted the spikes using old code).

path_name = '/local-data/DA-paper/fast-flicker/orn/cropped/';
[PID, ~, fA, paradigm, orn, fly, AllControlParadigms, ~, ~, spikes] = consolidateData(path_name,true);

% relabel all paradigms based on concentrations
conc = NaN*paradigm;
for i = 1:length(paradigm)
	a = strfind(AllControlParadigms(paradigm(i)).Name,'-');
	if ~isempty(a)
		z = strfind(AllControlParadigms(paradigm(i)).Name,'%');
		conc(i) = str2double(AllControlParadigms(paradigm(i)).Name(a+1:z-1));
	end
end

% shrink all data down to the essentials
PID = PID(35e3:55e3,~isnan(conc));
fA = fA(35e3:55e3,~isnan(conc));
orn = orn(~isnan(conc));
fly = fly(~isnan(conc));
spikes = spikes(~isnan(conc),35e4:55e4)';
conc = conc(~isnan(conc));

% make a color map
c = parula(9);
all_conc = unique(conc);

% make a time vector
time = 1e-3*(1:length(PID));

%% Stimulus Overview 
% In the following figure, we plot the stimulus, grouped by paradigm. 


figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
for i = 1:length(c)-1
	x = PID(:,conc==all_conc(i));
	if ~isvector(x)
		errorShade(time,mean(x,2),sem(x'),'Color',c(i,:));
	else
		plot(time,x,'Color',c(i,:));
	end
end

prettyFig()

if being_published
	snapnow
	delete(gcf)
end

%%
% Do all stimulus cases have the same variance? In the following plot we compare the stimulus distributions, ignoring for the moment the means. 

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
for i = 1:length(c) - 1
	stim = PID(:,conc==all_conc(i));

	if ~isvector(stim)
		for j = 1:width(stim)
			stim(:,j) = stim(:,j) - mean(stim(:,j));
		end
		[y,x] = hist(stim,50);
		y = y/sum(y);
		plot(x,mean(y,2),'Color',c(i,:))
	end
end

prettyFig()

if being_published
	snapnow
	delete(gcf)
end


%% Firing Gain
% In this section we back out filters on a trial-wise basis, and then plot the gain vs. the mean stimulus:

a = 1e3; z = length(PID);
[K,fp,gain,gain_err] = extractFilters(PID,fA,'use_cache',true,'a',a,'z',z);

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(mean(PID),gain,'k+')
x = mean(PID)'; 
ff = fit(x(~isnan(gain)),gain(~isnan(gain)),'power1');
plot(sort(x),ff(sort(x)),'r')
set(gca,'XScale','log','YScale','log')
xlabel('Mean Stimulus (V)')
ylabel('ORN Gain (Hz/V)')
prettyFig()

if being_published
	snapnow
	delete(gcf)
end

%%
% There seems to be a great deal of variability. Perhaps this is because we are including every single trial. What happens if we include only the first trial for each ORN? 

first_trials = true(length(conc),1);
for i = 2:length(conc)
	if conc(i) - conc(i-1) == 0  && orn(i)- orn(i-1) == 0
		first_trials(i) = false;
	end
end

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
x = mean(PID(:,first_trials));
y = gain(first_trials); 
x = x(~isnan(y)); y = y(~isnan(y));
plot(x,y,'k+')
ff = fit(x(:),y(:),'power1');
plot(sort(x),ff(sort(x)),'r')
set(gca,'XScale','log','YScale','log')
xlabel('Mean Stimulus (V)')
ylabel('ORN Gain (Hz/V)')
prettyFig()

if being_published
	snapnow
	delete(gcf)
end


%%
% What if we normalized the ORN gain to the ORN gain at the lowest dose? Perhaps this variability comes from ORN to ORN variability? 

lowest_gains = gain((conc == 3)' & first_trials);
norm_gain = gain;
for i = 1:max(orn)
	norm_gain(orn == i) = gain(orn==i)/lowest_gains(i);
end
norm_gain = norm_gain*nanmean(gain(conc==3));


figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(mean(PID(:,first_trials)),norm_gain(first_trials),'k+')
set(gca,'XScale','log','YScale','log')
xlabel('Mean Stimulus (V)')
ylabel('ORN Gain (Hz/V)')
prettyFig()

if being_published
	snapnow
	delete(gcf)
end

%%
% It's pretty clear that there's something very weird about the lowest dose. Let's skip that, and plot the rest. 

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
x = mean(PID(:,conc>3 & first_trials')); y = norm_gain(conc>3 & first_trials');
x = x(~isnan(y)); y = y(~isnan(y));
plot(x,y,'k+')
ff = fit(x(:),y(:),'power1');
plot(sort(x),ff(sort(x)),'r')
set(gca,'XScale','log','YScale','log')
xlabel('Mean Stimulus (V)')
ylabel('ORN Gain (Hz/V)')
prettyFig()

if being_published
	snapnow
	delete(gcf)
end

%% Kinetics
% In this section we look at the kinetics of the response. 

a = 1e3; z = length(PID);
tau = NaN(width(PID),1);
for i = 1:width(PID)
	x = PID(a:z,i);
	x = x - mean(x); x = x/std(x);
	y = fA(a:z,i);
	y = y - mean(y); y = y/std(y);
	temp = xcorr(x,y);
	[~,loc] = max(temp);
	tau(i) = z-a - loc;
end
tau(tau<0) = NaN;
tau(tau>1e3) = NaN;

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
for i = 2:max(all_conc)
	x = mean(mean(PID(:,conc == all_conc(i))));
	y = nanmean(tau(conc == all_conc(i)));
	ye = sem(tau(conc == all_conc(i)));
	errorbar(x,y,ye,'k')
end
title('PID\rightarrow Firing')
ylabel('Response delay (ms)')
xlabel('Mean Stimulus (V)')
set(gca,'YLim',[0 150])


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

