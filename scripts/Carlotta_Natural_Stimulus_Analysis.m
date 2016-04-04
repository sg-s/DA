% Carlotta_Natural_Stimulus_Analysis.m
% this script analyses naturalistic stimulus data from carlotta's paper, (fig 2)
% 
% created by Srinivas Gorur-Shandilya at 9:50 , 28 December 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

pHeader

%% Analysis of Carlotta's Naturalistic Stimulus
% In this document, we analyse data that Carlotta collected in Fig 2 of her paper.

% load the data

% find the longest trial
t_max = 0;
ntrials_1but = 0;

p = '/local-data/carlotta/2011_02_09_ab3A_evaporation/1but/*.mat';
allfiles = dir(p);
for i = 1:length(allfiles)
	clear PIDresp
	try
		load([fileparts(p),oss,allfiles(i).name])
		t_max = max([t_max length(PIDresp)]);
		ntrials_1but = ntrials_1but + width(PIDresp);
	catch
	end
end

% placeholders
PID_1but = zeros(t_max/10,ntrials_1but);

p = '/local-data/carlotta/2011_02_09_ab3A_evaporation/1but/*.mat';
allfiles = dir(p);
c = 1;
for i = 1:length(allfiles)
	clear PIDresp
	try
		load([fileparts(p),oss,allfiles(i).name]);
		temp = PIDresp(1:10:end,:);
		PID_1but(1:length(temp),c:c+width(temp)-1) = temp;
		c = c + width(temp) + 1;
	catch
	end
end

% throw out traces that are less than 20s long
rm_this = (PID_1but(100e3,:) == 0);
PID_1but(:,rm_this) = [];
PID_1but(PID_1but==0) = NaN;

% throw out small, short traces
rm_this = max(PID_1but) < .2;
PID_1but(:,rm_this) = [];

% remove baseline from first 5 seconds
for i = 1:width(PID_1but)
	PID_1but(:,i) = PID_1but(:,i) - nanmean(PID_1but(1:5e3,i));
end

time = 1e-3*(1:length(PID_1but));

%%
% This is what the raw data looks like. Each trace is a different trial. 

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
plot(time,PID_1but)
xlabel('Time (s)')
ylabel('Stimulus (V)')
title('Methyl butyrate')
prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%%
% We now plot the stimulus distributions for two time windows along the stimulus. 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
for i = 1:width(PID_1but)
	X = PID_1but(20e3:25e3,i);
	hx = linspace(0,0.3,100);
	hy = histcounts(X,hx);
	hy = hy/sum(hy);
	hx = hx(1:end-1) + mean(diff(hx));
	plot(hx,hy)
end
title('20-25s')
xlabel('Stimulus (V)')
ylabel('Probability')

subplot(1,2,2), hold on
for i = 1:width(PID_1but)
	X = PID_1but(80e3:85e3,i);
	hx = linspace(0,0.3,100);
	hy = histcounts(X,hx);
	hy = hy/sum(hy);
	hx = hx(1:end-1) + mean(diff(hx));
	plot(hx,hy)
end
title('80-85s')
xlabel('Stimulus (V)')
prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%%
% It looks like the stimulus gets smaller and narrower as time progresses. Now, we look at the correlation between the mean and the variance in these stimuli over various timescales.

figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[800 800]); hold on
for i = 4:-1:1
	ax(i) = subplot(2,2,i); hold on
	xlabel('\mu_{Stimulus} (V)')
	ylabel('\sigma_{Stimulus} (V)')
	set(gca,'XLim',[0 0.3],'YLim',[0 0.3])
	plot([0 0.3],[0 0.3],'k--')
end

% now show that the variance and the mean are correlated 
all_block_sizes = [10 500 1e3 5e3];
clear l r2
r2 = NaN*all_block_sizes;

for j = 1:width(PID_1but)
	for i = 1:length(all_block_sizes)
		temp = nonnans(PID_1but(:,j));
		temp(temp<0) = [];
		z = floor(length(temp)/all_block_sizes(i))*all_block_sizes(i);
		temp = temp(1:z);
		temp = reshape(temp,all_block_sizes(i),length(temp)/all_block_sizes(i));
		plot(ax(i),mean(temp),std(temp),'Marker','.','Color',[.5 .5 .5],'LineStyle','none')
		r2(i) = rsquare(mean(temp),std(temp));
		title(ax(i),[oval(all_block_sizes(i)) 'ms, r^2 = ' oval(r2(i))])
	end
end

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end



%% Version Info
%
pFooter;


