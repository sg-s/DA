% Analysis of OBP28a Knockout
% 
% created by Srinivas Gorur-Shandilya at 10:52 , 20 August 2015. Contact me at http://srinivas.gs/contact/
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


%% Analysis of OBP28a Knockout
% In this document, we analyse the responses of ab3A in a OBP28a knockout (using Crispr). 

%% LFP Filter Analysis
% First, we analyse the LFP. In the following figure, we extract LN models for the LFP in all cases, and compare the OBP knockout to a control genotype. The odour used is ethyl acetate. 

[PID, LFP, fA, paradigm, orn] = consolidateData('/local-data/obp/ko',1);
geno = ones(length(orn),1);

[PID2, LFP2, fA2, paradigm2, orn2] = consolidateData('/local-data/obp/wcs',1);
geno2 = 2*ones(length(orn2),1);

PID = [PID PID2]; clear PID2
LFP = [LFP LFP2]; clear LFP2
fA = [fA fA2]; clear fA2
paradigm = [paradigm paradigm2]; clear paradigm2
orn = [orn orn2]; clear orn2
geno = [geno; geno2]; clear geno2
paradigm = paradigm(:);
orn = orn(:);

% filter the LFP
filteredLFP = LFP;
for i = 1:width(LFP)
	if paradigm(i)==2
		a = find(~isnan(LFP(:,i)),1,'first');
		z = find(~isnan(LFP(:,i)),1,'last');
		if isempty(a)
			a = 1;
		end
		if isempty(z)
			z = length(LFP);
		end
		filteredLFP(a:z,i) = 10*filter_trace(LFP(a:z,i),1000,10); % now in mV
	end
end

% K -- PID -> LFP filter
K = cache(DataHash([PID; LFP]));
if isempty(K);
	K = NaN(1e3,length((orn)));
	for i = 1:length((orn))
		textbar(i,length(orn))
		if paradigm(i) == 2
			resp = (filteredLFP(10e3:end,i));
			stim = (PID(10e3:end,i));
			stim(1:400) = [];
			resp(end-399:end) = [];
			stim = stim - mean(stim);
			stim =  stim/std(stim);
			resp  =resp - mean(resp);
			resp = resp/std(resp);
			temp = fitFilter2Data(stim,resp,'reg',1,'filter_length',1399);
			% throw out 200ms on either end
			temp(1:200) = [];
			temp(end-199:end) = [];
			K(:,i) = temp;
		end
	end
end
% cache
cache(DataHash([PID; filteredLFP]),K);

% make linear predictions everywhere
time = 1e-3*(1:length(filteredLFP));
LFP_pred = NaN*LFP;
for i = 1:length((orn))
	if paradigm(i) == 2
		filtertime = 1e-3*(1:1e3)-.2;
		LFP_pred(:,i) = convolve(time,PID(:,i),K(:,i),filtertime);
	end
end

% fit nonlinearities everywhere
clear ff
for i = 1:length((orn))
	if paradigm(i) == 2
		x = LFP_pred(10e3:end-1e3,i);
		y = filteredLFP(10e3:end-1e3,i); y = y-mean(y);
		ff(i).LFP = fit(x(:),y(:),'poly5');
	end
end

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
clear l 
l(1) = errorShade(filtertime,mean2(K(:,geno==1)),sem(K(:,geno==1)),'Color',[1 0 0]);
l(2) = errorShade(filtertime,mean2(K(:,geno==2)),sem(K(:,geno==2)),'Color',[0 0 1]);
legend(l,{'CRISPR KO','wCS'})
xlabel('Lag (s)')
ylabel('Filter Amplitude')


subplot(1,2,2), hold on
x = -.8:0.005:.4;
y = NaN(length(x),length(orn));
for i = 1:length((orn))
	if paradigm(i) == 2
		y(:,i) = ff(i).LFP(x(:));
	end
end
clear l 
l(1) = errorShade(x,mean2(y(:,geno==1)),sem(y(:,geno==1)),'Color',[1 0 0]);
l(2) = errorShade(x,mean2(y(:,geno==2)),sem(y(:,geno==2)),'Color',[0 0 1]);
xlabel('Filter Output')
ylabel('LFP')

PrettyFig()

if being_published
	snapnow
	delete(gcf)
end


%% LFP Pulse Analysis
% We now look for differences in the LFP to a pulse of odour. 

% remove baseline
for i = 1:length(orn)
	if paradigm(i)==3
		filteredLFP(:,i) = filteredLFP(:,i) - mean(filteredLFP(1:5e3,i));
	end
end

figure('outerposition',[0 0 700 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
clear l 
l(1) = errorShade(time(1:end-10),mean2(filteredLFP(1:end-10,geno==1 & paradigm == 3)),sem(filteredLFP(1:end-10,geno==1 & paradigm == 3)),'Color',[1 0 0]);
l(2) = errorShade(time(1:end-10),mean2(filteredLFP(1:end-10,geno==2 & paradigm == 3)),sem(filteredLFP(1:end-10,geno==2 & paradigm == 3)),'Color',[0 0 1]);
legend(l,{'CRISPR KO','wCS'},'location','SouthEast')
xlabel('Lag (s)')
ylabel('LFP (mV)')

PrettyFig()

if being_published
	snapnow
	delete(gcf)
end

%% Firing Rate Pulse Analysis
% We now look at the firing rates in response to a long pulse. 

% remove baseline
for i = 1:length(orn)
	if paradigm(i)==3
		fA(:,i) = fA(:,i) - mean(fA(1:5e3,i));
	end
end

figure('outerposition',[0 0 700 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
clear l 
l(1) = errorShade(time(1:end-10),mean2(fA(1:end-10,geno==1 & paradigm == 3)),sem(fA(1:end-10,geno==1 & paradigm == 3)),'Color',[1 0 0]);
l(2) = errorShade(time(1:end-10),mean2(fA(1:end-10,geno==2 & paradigm == 3)),sem(fA(1:end-10,geno==2 & paradigm == 3)),'Color',[0 0 1]);
legend(l,{'CRISPR KO','wCS'},'location','NorthEast')
xlabel('Lag (s)')
ylabel('Firing Rate (Hz)')

PrettyFig()

if being_published
	snapnow
	delete(gcf)
end

%%
% In the following figure, we rescale each of the firing rates by the peak, and then plot it, to compare the dynamical nature of the firing rate. 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,3,1:2), hold on
clear l 
temp = mean2(fA(1:end-10,geno==1 & paradigm == 3));
temp = temp/max(temp);
l(1) = plot(time(1:end-10),temp,'Color',[1 0 0]);
temp = mean2(fA(1:end-10,geno==2 & paradigm == 3));
temp = temp/max(temp);
l(2) = plot(time(1:end-10),temp,'Color',[0 0 1]);
legend(l,{'CRISPR KO','wCS'},'location','NorthEast')
xlabel('Time (s)')
ylabel('Firing Rate (norm)')

subplot(1,3,3), hold on
temp = mean2(fA(1:end-10,geno==1 & paradigm == 3));
temp = temp/max(temp);
l(1) = plot(time(1:end-10),temp,'Color',[1 0 0]);
temp = mean2(fA(1:end-10,geno==2 & paradigm == 3));
temp = temp/max(temp);
l(2) = plot(time(1:end-10),temp,'Color',[0 0 1]);
xlabel('Time (s)')
ylabel('Firing Rate (norm)')
set(gca,'XLim',[4 9])

PrettyFig()

if being_published
	snapnow
	delete(gcf)
end

%% Firing Rate Filter Analysis
% We now extract filters for the firing rate of each neuron, and compare the two genotypes. 



% K2 -- PID -> fA filter
K2 = cache(DataHash([PID; fA]));
if isempty(K2);
	K2 = NaN(1e3,length((orn)));
	for i = 1:length((orn))
		textbar(i,length(orn))
		if paradigm(i) == 2
			resp = (fA(10e3:end,i));
			stim = (PID(10e3:end,i));
			stim(1:400) = [];
			resp(end-399:end) = [];
			stim = stim - mean(stim);
			stim =  stim/std(stim);
			resp  =resp - mean(resp);
			resp = resp/std(resp);
			temp = fitFilter2Data(stim,resp,'reg',1,'filter_length',1399);
			% throw out 200ms on either end
			temp(1:200) = [];
			temp(end-199:end) = [];
			K2(:,i) = temp;
		end
	end
end
% cache
cache(DataHash([PID; fA]),K2);

% make linear predictions everywhere
time = 1e-3*(1:length(fA));
fp = NaN*LFP;
for i = 1:length((orn))
	if paradigm(i) == 2
		filtertime = 1e-3*(1:1e3)-.2;
		fp(:,i) = convolve(time,PID(:,i),K2(:,i),filtertime);
	end
end

% fit nonlinearities everywhere--this is slow, so cached
ff = cache(DataHash([fp fA]));
if isempty(ff)
	for i = 1:length((orn))
		textbar(i,length(orn))
		if paradigm(i) == 2
			x = fp(10e3:end-1e3,i);
			y = fA(10e3:end-1e3,i);

			ft = fittype( 'hill_fit(x,A,k,n,offset)' );
			ff(i).fA = fit(x(:),y(:),ft,'StartPoint',[100, 50, 2,3],'Lower',[1 eps 1 -10],'Upper',[1e4 1e3 10 10]);
		end
	end
end
cache(DataHash([fp fA]),ff);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
clear l 
l(1) = errorShade(filtertime,mean2(K2(:,geno==1)),sem(K2(:,geno==1)),'Color',[1 0 0]);
l(2) = errorShade(filtertime,mean2(K2(:,geno==2)),sem(K2(:,geno==2)),'Color',[0 0 1]);
legend(l,{'CRISPR KO','wCS'})
xlabel('Lag (s)')
ylabel('Filter Amplitude')


subplot(1,2,2), hold on
x = -.5:0.005:1;
y = NaN(length(x),length(orn));
for i = 1:length((orn))
	if paradigm(i) == 2
		y(:,i) = ff(i).fA(x(:));
	end
end
clear l 
l(1) = errorShade(x,mean2(y(:,geno==1)),sem(y(:,geno==1)),'Color',[1 0 0]);
l(2) = errorShade(x,mean2(y(:,geno==2)),sem(y(:,geno==2)),'Color',[0 0 1]);
xlabel('Filter Output')
ylabel('Firing Rate (Hz)')

PrettyFig()

if being_published
	snapnow
	delete(gcf)
end

%% Fast Gain Control Analysis
% We now check for fast gain control in the two genotypes, and see if there are any differences between the two. 

% make linear predictions everywhere
time = 1e-3*(1:length(fA)); time = time(:);
fp = NaN*LFP;
for i = 1:length((orn))
	if paradigm(i) == 2
		filtertime = 1e-3*(1:1e3)-.2;
		fp(:,i) = convolve(time,PID(:,i),K2(:,i),filtertime);

		% correct for some trivial scaling
		rm_this = isnan(fp(:,i)) | isnan(fA(:,i)) | time < 10;
		x = fp(:,i); y = fA(:,i);
		x(rm_this) = []; y(rm_this) = [];
		temp = fit(x,y,'poly1');
		fp(:,i) = temp(fp(:,i));
	end
end


% gain analysis -- Linear model
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on

% first we do the firing rates
clear ph ax
ax(1) = subplot(1,2,1); hold on
ax(2) = subplot(1,2,2); hold on

history_lengths = logspace(-1,1,30);

for g = 1:2
	for i = 1:max(orn(geno==g))

		ph(4) = ax(g);

		resp = mean2(fA (:,orn==i & geno == g & paradigm == 2));
		stim = mean2(PID(:,orn==i & geno == g & paradigm == 2));
		pred = mean2(fp (:,orn==i & geno == g & paradigm == 2));

		try
			stim = stim(20e3:55e3);
			pred = pred(20e3:55e3);
			resp = resp(20e3:55e3);

			[p,~,~,~,~,history_lengths]=GainAnalysisWrapper2('response',resp,'prediction',pred,'stimulus',stim,'time',1e-3*(1:length(resp)),'ph',ph,'history_lengths',history_lengths,'use_cache',1,'engine',@GainAnalysis5);
		end
	end
end


for g = 1:2
	h=get(ax(g),'Children');
	rm_this = [];
	for i = 1:length(h)
		try
			if  strcmp(get(h(i),'LineStyle'),'-.')
				rm_this = [rm_this i];
			end
		catch
		end
	end
	delete(h(rm_this))

	set(ax(g),'XLim',[.1 10],'YLim',[.4 2.5])
end

title(ax(1),'CRISPR KO')
title(ax(2),'wCS')

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
	unix(strjoin({'tag -a published',which(mfilename)}));
end
