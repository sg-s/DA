% PTX.m
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


%% Data Overview
% In this section we compare the stimulus and LFP in RNAi-knockdown flies with parental controls. 

[alldata(1).PID, alldata(1).LFP, alldata(1).fA, alldata(1).paradigm, alldata(1).orn, alldata(1).fly, alldata(1).AllControlParadigms, alldata(1).paradigm_hashes] = consolidateData('/local-data/DA-paper/ptx/uas-gal4/',false);
alldata(1).genotype = 'UAS-PTX-22a-gal4';

[alldata(2).PID, alldata(2).LFP, alldata(2).fA, alldata(2).paradigm, alldata(2).orn, alldata(2).fly, alldata(2).AllControlParadigms, alldata(2).paradigm_hashes] = consolidateData('/local-data/DA-paper/ptx/uas-control/',false);
alldata(2).genotype = 'UAS-control';

[alldata(3).PID, alldata(3).LFP, alldata(3).fA, alldata(3).paradigm, alldata(3).orn, alldata(3).fly, alldata(3).AllControlParadigms, alldata(3).paradigm_hashes] = consolidateData('/local-data/DA-paper/g-alpha/rnai/gal4-control/',false);
alldata(3).genotype = 'GAL4-control';


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
for ai = 1:length(alldata)
	subplot(2,3,ai), hold on
	ss = 10;
	for i = 1:length(all_paradigm_hashes)
		this_hash = all_paradigm_hashes{i};
		this_paradigm = find(strcmp(alldata(ai).paradigm_hashes,this_hash));
		plot_this = alldata(ai).PID(1:ss:end,alldata(ai).paradigm==this_paradigm);
		time = 1e-2*(1:length(plot_this));
		try
			plot(time,mean(plot_this,2),'Color',c(this_paradigm,:),'LineWidth',2);
		catch
		end
	end
	ylabel('Stimulus (V)')
	xlabel('Time (s)')
	set(gca,'YLim',[0 2.5])
	title(alldata(ai).genotype)

	subplot(2,3,ai+length(alldata)), hold on
	for i = 1:length(all_paradigm_hashes)
		this_hash = all_paradigm_hashes{i};
		this_paradigm = find(strcmp(alldata(ai).paradigm_hashes,this_hash));
		plot_this = alldata(ai).LFP(1:ss:end,alldata(ai).paradigm==this_paradigm);
		time = 1e-2*(1:length(plot_this));
		try
			plot(time,mean(plot_this,2),'Color',c(this_paradigm,:),'LineWidth',2);
		catch
		end
	end
	ylabel('LFP (mV)')
	xlabel('Time (s)')
	set(gca,'YLim',[-4 .5])
end

prettyFig('fs=14;');

if being_published
	snapnow
	delete(gcf)
end


%% LFP Gain
% In this section we analyse the gain in the LFP and compare it to the wild type. 


a = 20e3; z = 50e3;
for ai = 1:length(alldata)
	[alldata(ai).K1,alldata(ai).LFP_pred,alldata(ai).LFP_gain,alldata(ai).LFP_gain_err] = extractFilters(alldata(ai).PID,alldata(ai).filtered_LFP,'use_cache',true,'a',a,'z',z);
end


figure('outerposition',[0 0 600 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
ai = 2;
clear l
for i = 1:width(alldata(ai).LFP)
	l(ai) = plot(mean(alldata(ai).PID(a:z,i)),alldata(ai).LFP_gain(i),'+','Color','g');
end
ai = 3;
for i = 1:width(alldata(ai).LFP)
	l(ai) = plot(mean(alldata(ai).PID(a:z,i)),alldata(ai).LFP_gain(i),'+','Color','k');
end
ai = 1;
for i = 1:width(alldata(ai).LFP)
	l(ai) = plot(mean(alldata(ai).PID(a:z,i)),alldata(ai).LFP_gain(i),'+','Color','r');
end

set(gca,'XScale','log','YScale','log')
xlabel('Mean Stimulus (V)')
ylabel('LFP Gain (mV/V)')
legend(l,{alldata.genotype})

set(gca,'XLim',[.1 10],'YLim',[.05 10])

prettyFig;

if being_published
	snapnow
	delete(gcf)
end

%% Kinetics 
% In this section, we check if there are any kinetics changes in the mutant, and compare them to the controls. 

% a = 30e3; z = 50e3;
% for ai = 1:length(alldata)
% 	alldata(ai).tau = NaN(width(alldata(ai).PID),1);
% 	for i = 1:width(alldata(ai).PID)
% 		x = alldata(ai).PID(a:z,i);
% 		x = x - mean(x); x = x/std(x);
% 		y = alldata(ai).filtered_LFP(a:z,i);
% 		y = y - mean(y); y = y/std(y);
% 		y = -y;
% 		temp = xcorr(x,y);
% 		[~,loc] = max(temp);
% 		alldata(ai).tau(i) = 20e3 - loc;
% 	end
% end

% figure('outerposition',[0 0 600 500],'PaperUnits','points','PaperSize',[600 500]); hold on
% ai = 1;
% for i = 1:max(alldata(ai).paradigm)
% 	mean_stim = mean(mean(alldata(ai).PID(a:z,alldata(ai).paradigm==i)));
% 	mean_tau = mean(alldata(ai).tau(alldata(ai).paradigm==i));
% 	tau_err = sem(alldata(ai).tau(alldata(ai).paradigm==i));
% 	errorbar(mean_stim,mean_tau,tau_err,'r')
% end

% ai = 2;
% for i = 1:max(alldata(ai).paradigm)
% 	mean_stim = mean(mean(alldata(ai).PID(a:z,alldata(ai).paradigm==i)));
% 	mean_tau = mean(alldata(ai).tau(alldata(ai).paradigm==i));
% 	tau_err = sem(alldata(ai).tau(alldata(ai).paradigm==i));
% 	errorbar(mean_stim,mean_tau,tau_err,'g')
% end

% ai = 3;
% for i = 1:max(alldata(ai).paradigm)
% 	mean_stim = mean(mean(alldata(ai).PID(a:z,alldata(ai).paradigm==i)));
% 	mean_tau = mean(alldata(ai).tau(alldata(ai).paradigm==i));
% 	tau_err = sem(alldata(ai).tau(alldata(ai).paradigm==i));
% 	errorbar(mean_stim,mean_tau,tau_err,'k')
% end


% xlabel('Mean Stimulus (V)')
% ylabel('Response delay (ms)')

% prettyFig;

% if being_published
% 	snapnow
% 	delete(gcf)
% end

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

%%
% This file has the following external dependencies:
showDependencyHash(mfilename);

t = toc;

%% 
% This document was built in: 
disp(strcat(oval(t,3),' seconds.'))

% tag the file as being published 

if being_published
	unix(['tag -a published ',which(mfilename)]);
	unix(['tag -r publish-failed ',which(mfilename)]);
end

