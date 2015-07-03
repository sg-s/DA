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

%% Stimulus Characterisation
% First, we show that we are able to deliver Gaussian-distributed odour stimuli, and that we are able to vary the means of the distributions of these stimuli. 

PID = [];
paradigm = [];
fly = [];
AllControlParadigms = struct;
AllControlParadigms.Name = '';
AllControlParadigms.Outputs = [];
AllControlParadigms(1) = [];
paradigm_hashes = {};

allfiles = dir('/local-data/DA-paper/LFP-MSG/*.mat');
for i = 1:length(allfiles)
	load(strcat('/local-data/DA-paper/LFP-MSG/',allfiles(i).name));
	for j = 1:length(data)
		if ~isempty(data(j).PID)
			clear this_paradigm 
			% figure out which control paradigm this is
			this_hash = DataHash(ControlParadigm(j));
			if isempty(find(strcmp(this_hash,paradigm_hashes)))
				AllControlParadigms(end+1) = ControlParadigm(j);
				this_paradigm = length(AllControlParadigms);
				paradigm_hashes{end+1} = this_hash;
			else
				this_paradigm = find(strcmp(this_hash,paradigm_hashes));
			end


			this_PID = data(j).PID;
			if length(spikes) < j
			else
				rm_this = [];
				try
					rm_this = find(spikes(j).discard);
				end
				if ~isempty(rm_this)
					this_PID(rm_this,:) = [];
				else

				end

				this_PID = this_PID(:,1:10:end)';
				PID = [PID this_PID];
				paradigm = [paradigm  this_paradigm*ones(1,width(this_PID))];
				fly = [fly  i*ones(1,width(this_PID))];
			end
		end
	end
end

return

% sort the paradigms sensibily
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



%%
% The following figure shows the distribution of the inputs, for the terminal 40seconds of each 60second presentation. 
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
c = parula(length(unique(paradigm)));
for i = 1:length(c)
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

