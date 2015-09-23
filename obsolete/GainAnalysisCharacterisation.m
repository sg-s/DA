% GainAnalysisCharacterisation.m
% careful characterisation of the gain analysis method
% created by Srinivas Gorur-Shandilya at 10:22 , 19 July 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

% this code determines if this function is being called by publish() or not
calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end
tic

%% Negative Controls
% In this section, we check that our gain analysis methods never suggest that fast gain adaptation is happening in synthetic data constructed from LN models (which, by definition, should not show fast gain adaptation). 

%%
% The following figure shows the models we use to generate the synthetic data. The models are made from a fixed bilobed filter and a range of saturating non-linearities.  

K = (filter_alpha2(50,100,1,.4,1:1000));
A = 100;
s = fliplr([25 50 100 1000]);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
plot(1:1000,K)
xlabel('Filter lag (ms)')
title('Filter')
subplot(1,2,2), hold on
c = parula(length(s)+1);
for i = 1:length(s)
	x = [0:0.01:.99 1:.1:100];
	y = hill([A s(i) 1],x);
	y = y/max(y);
	plot(x,y,'Color',c(i,:))
end
title('Output Nonlinearity')
xlabel('Filter Output (a.u.)')

PrettyFig()

if being_published
	snapnow
	delete(gcf)
end

%%
% Now, we use these models to generate synthetic data from the stimulus we use in the actual experiments:

load('/local-data/DA-paper/large-variance-flicker/LFP/consolidated_data.mat')
PID = mean2(cached_data.PID) + .2;
clear cached_data
time = 1e-3*(1:length(PID));

figure('outerposition',[0 0 700 700],'PaperUnits','points','PaperSize',[700 700]); hold on
subplot(2,1,1), hold on
plot(time,PID,'k')
xlabel('Time (s)')
ylabel('Stimulus (V)')

noise_scale = .1;
subplot(2,1,2), hold on
resp = NaN(length(PID),length(s));
for i = 1:length(s)
	rs = RandStream('mt19937ar','Seed',1);
	RandStream.setGlobalStream(rs);
	resp(:,i) = filter(K,1,noise_scale*randn(length(PID),1)+PID);
end


for i = 1:length(s)
	resp(:,i) = hill([A s(i) 1],resp(:,i));
	resp(:,i) = (100*resp(:,i))/hill([A s(i) 1],100);
	plot(time,resp(:,i),'Color',c(i,:))
end
xlabel('Time (s)')
ylabel('Response (Hz)')

PrettyFig()

if being_published
	snapnow
	delete(gcf)
end

%%
% We now perform gain analysis on each of these synthetic data sets, and check if our analysis detects a signature of what we think is fast gain control. In this analysis, we assume that we know the filter (or our filter estimation technique is perfect). 

% first, make linear predictions for everything
pred = NaN*resp;
for i = 1:length(s)
	pred(:,i) = filter(K,1,PID);
	% correct for some trivial scaling
	ff = fit(pred(:,i),resp(:,i),'poly1');
	pred(:,i) =  ff(pred(:,i));
end

figure('outerposition',[0 0 1400 800],'PaperUnits','points','PaperSize',[1400 800]); hold on
for i = 1:length(s)
	clear ph
	ph(3) = subplot(2,length(s),i); hold on
	ph(4) = subplot(2,length(s),length(s)+i); hold on

	hl_min = .1;
	hl_max = 10;
	history_lengths = logspace(log10(hl_min),log10(hl_max),30);

	history_lengths = findValidHistoryLengths(1e-3,PID,pred(:,i),resp(:,i),30,.33);

	[p,~,~,~,~,history_lengths]=GainAnalysisWrapper2('response',resp(:,i),'prediction',pred(:,i),'stimulus',PID,'time',1e-3*(1:length(resp)),'ph',ph,'history_lengths',history_lengths,'use_cache',1,'engine',@GainAnalysis5,'example_history_length',history_lengths(11));

end

PrettyFig()

if being_published
	snapnow
	delete(gcf)
end

%%
% There seems to be some false positive rate, but it's not too bad. 

%%
% What if we first back out the filters from the data, and use those in our gain analysis? The following figure shows the backed-out filters for the dataset:

Khat = NaN(1e3,length(s));

for i = 1:length(s)
	x = PID; y = resp(:,i);
	x(1:400) = [];
	y(end-399:end) = [];
	this_K = FitFilter2Data(x,y,[],'filter_length=1399;','reg=1;');
	Khat(:,i) = this_K(401:end);
end

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
c = parula(length(s)+1);
for i = 1:length(s)
	plot(1:1000,Khat(:,i),'Color',c(i,:))
end
xlabel('Filter Lag (ms)')
title('Reconstructed Effective Filters')

PrettyFig()

if being_published
	snapnow
	delete(gcf)
end

%%
% We now use these backed out filters to make linear predictions of the data, and then run gain analysis on these predictions: 

pred = NaN*resp;
for i = 1:length(s)
	pred(:,i) = filter(Khat(:,i),1,PID);
	% correct for some trivial scaling
	ff = fit(pred(:,i),resp(:,i),'poly1');
	pred(:,i) =  ff(pred(:,i));
end

figure('outerposition',[0 0 1400 800],'PaperUnits','points','PaperSize',[1400 800]); hold on
for i = 1:length(s)
	clear ph
	ph(3) = subplot(2,length(s),i); hold on
	ph(4) = subplot(2,length(s),length(s)+i); hold on

	hl_min = .1;
	hl_max = 10;
	history_lengths = logspace(log10(hl_min),log10(hl_max),30);

	history_lengths = findValidHistoryLengths(1e-3,PID,pred(:,i),resp(:,i),30,.33);

	[p,~,~,~,~,history_lengths]=GainAnalysisWrapper2('response',resp(:,i),'prediction',pred(:,i),'stimulus',PID,'time',1e-3*(1:length(resp)),'ph',ph,'history_lengths',history_lengths,'use_cache',1,'engine',@GainAnalysis5,'example_history_length',history_lengths(11));

end

PrettyFig()

if being_published
	snapnow
	delete(gcf)
end


%% Positive Controls
% Now, we check the ability of the analysis method to detect fast gain control in data from a model that explicitly assumes fast gain control (i.e., the DA model). The following figure shows the parameters of the DA model used to generate synthetic data. We also vary the B parameter (which has no visible effects on the filters, but varies the contribution of the $K_z$ filter)

p = getModelParameters('DAmodelv2');
p(2:4) = p;
A = [1 5 10 20];
B = [0 5 10 20];
for i = 1:length(p)
	p(i).A = A(i);
	p(i).B = B(i);
	p(i).C = .3;
	p(i).s0 = 0;
	p(i).tau_y = 50;
	p(i).tau_z = 100;
end

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
[Ky,Kz] = DA_Filters(p(1),1:1000);
plot(1:1000,Ky,'b')
plot(1:1000,Kz,'r')
xlabel('Filter Lag (ms)')
ylabel('Filter Amplitude')
legend({'K_y','K_z'})
PrettyFig()

if being_published
	snapnow
	delete(gcf)
end

%%
% We now look at the response generated from this model:

figure('outerposition',[0 0 700 700],'PaperUnits','points','PaperSize',[700 700]); hold on
subplot(2,1,1), hold on
plot(time,PID,'k')
xlabel('Time (s)')
ylabel('Stimulus (V)')

noise_scale = .1;
subplot(2,1,2), hold on
resp = NaN(length(PID),length(s));
c = parula(length(B)+1);
for i = 1:length(B)
	rs = RandStream('mt19937ar','Seed',1);
	RandStream.setGlobalStream(rs);
	resp(:,i) = 100*DAModelv2(noise_scale*randn(length(PID),1)+PID,p(i));
	plot(time,resp(:,i),'Color',c(i,:))
end

xlabel('Time (s)')
ylabel('Response (Hz)')

PrettyFig()

if being_published
	snapnow
	delete(gcf)
end

%%
% Now, we back out linear filters from this synthetic data:

Khat = NaN(1e3,length(B));

for i = 1:length(B)
	Khat(:,i) = FitFilter2Data(PID,resp(:,i),[],'filter_length=999;','reg=1;');
end

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
c = parula(length(B)+1);
for i = 1:length(B)
	plot(1:1000,Khat(:,i),'Color',c(i,:))
end
xlabel('Filter Lag (ms)')
title('Reconstructed Effective Filters')

PrettyFig()

if being_published
	snapnow
	delete(gcf)
end

%%
% We can see that as we turn up the $\beta$ parameter, the effective filter turns from single-lobed to a bi-lobed differentiating filter. Now, we check if we see evidence for fast gain control using our gain analysis methods:

% first, make linear predictions for everything

pred = NaN(length(PID),length(B));
for i = 1:length(B)
	pred(:,i) = filter(Khat(:,i),1,PID);
	% correct for some trivial scaling
	ff = fit(pred(:,i),resp(:,i),'poly1');
	pred(:,i) =  ff(pred(:,i));
end

figure('outerposition',[0 0 1400 800],'PaperUnits','points','PaperSize',[1400 800]); hold on
for i = 1:length(B)
	clear ph
	ph(3) = subplot(2,length(B),i); hold on
	ph(4) = subplot(2,length(B),length(B)+i); hold on

	hl_min = .1;
	hl_max = 10;
	history_lengths = logspace(log10(hl_min),log10(hl_max),30);

	history_lengths = findValidHistoryLengths(1e-3,PID,pred(:,i),resp(:,i),30,.33);

	GainAnalysisWrapper2('response',resp(1e3:end-1e3,i),'prediction',pred(1e3:end-1e3,i),'stimulus',PID(1e3:end-1e3),'time',1e-3*(1:length(PID(1e3:end-1e3))),'ph',ph,'history_lengths',history_lengths,'use_cache',1,'engine',@GainAnalysis5,'example_history_length',history_lengths(11));
	title(ph(3),['\beta=',oval(p(i).B)])

end

PrettyFig()

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

% tag the file as being published 
% add homebrew path
path1 = getenv('PATH');
path1 = [path1 ':/usr/local/bin'];
setenv('PATH', path1);

if being_published
	unix(strjoin({'tag -a published ',which(mfilename)}));
end
