% LVF_LFP.m
% recreates fig3 with LFP recordings
% 
% created by Srinivas Gorur-Shandilya at 11:52 , 14 July 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


pHeader;

%% 
% In this document, we attempt to recreate the results from Fig 3 showing fast gain changes with a widely distributed stimulus. 

%% Stimulus 
% First, we check that the stimulus is the same as the one we used six months ago (in LargeVarianceFlickering.m). The following plot shows the old stimulus (blue) and the new stimulus (red). 

% load old stimulus
load(getPath(dataManager,'b5169f03f72690c1ac6918632fb1f998'))
PID = data(4).PID;

load(getPath(dataManager,'baccc203cb85f35cec96172932275a26'))
PID = vertcat(PID,data(4).PID);
PID = PID';
PID = PID(1:10:end,:);
% remove baseline
baseline = mean(PID(1:300,1));
PID = PID - baseline;

% load new stimulus
load(getPath(dataManager,'d779046a10c49b570aba6d8bb665ea62'))
newPID = data(4).PID';
newPID = newPID(1:10:end,:);
time = 1e-3*(1:length(PID));

% remove baseline
baseline = mean(newPID(1:300,1));
newPID = newPID - baseline;

% plot
clear l
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
l(1) = errorShade(time(1:10:end),mean(PID(1:10:end,:),2),sem(PID(1:10:end,:)'),'Color',[0 0 1]);
l(2) = errorShade(time(1:10:end),mean(newPID(1:10:end,:),2),sem(newPID(1:10:end,:)'),'Color',[1 0 0]);
xlabel('Time (s)')
ylabel('PID (V)')
set(gca,'XLim',[20 60])
legend(l,{'January 28 2015','July 14 2015'})

prettyFig;

if being_published
	snapnow
	delete(gcf)
end

%%
% It looks very similar, but smaller, which is expected, as the PID gets less sensitive over time. What if we rescale it? 

ff = fit(mean(newPID,2),mean(PID,2),'poly1');

clear l
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
l(1) = errorShade(time(1:10:end),mean(PID(1:10:end,:),2),sem(PID(1:10:end,:)'),'Color',[0 0 1]);
l(2) = errorShade(time(1:10:end),ff(mean(newPID(1:10:end,:),2)),sem(newPID(1:10:end,:)'),'Color',[1 0 0]);
xlabel('Time (s)')
ylabel('PID (V)')
legend(l,{'January 28 2015','July 14 2015 Rescaled'})
set(gca,'XLim',[20 60])
prettyFig;

if being_published
	snapnow
	delete(gcf)
end

%% Summary of Data
% How does the neuron respond to this stimulus? Does the LFP and/or firing rates show evidence for fast gain control? In the following figure, we plot the stimulus, LFP and the response from all the neurons in the dataset. The shading shows the standard error of the mean. Each neuron is in a different colour. The stimulus, the LFP and the firing rates are all highly reproducible. 

p = getPath(dataManager,'c8dc5353a75ce5fcdcfa13139c716bd8');
[PID, LFP, fA, paradigm, orn, AllControlParadigms, paradigm_hashes,sequence] = consolidateData(p,1);

% set to NaN firing rates that are 0
fA(:,max(fA) == 0) = NaN;


% remove baseline
PID = PID - min(min(PID));

% throw out trials where we didn't record the LFP, for whatever reason
not_LFP = 0*orn;
for i = 1:width(LFP)
	not_LFP(i) = abs(mean(LFP(:,i)));
end
LFP(:,not_LFP< 0.5) = NaN;

% filter the LFP
filteredLFP = LFP;
for i = 1:width(LFP)
	a = find(~isnan(LFP(:,i)),1,'first');
	z = find(~isnan(LFP(:,i)),1,'last');
	if isempty(a)
		a = 1;
	end
	if isempty(z)
		z = length(LFP);
	end
	LFP(:,i) = LFP(:,i) - nanmean(LFP(:,i));
	try
		filteredLFP(a:z,i) = 10*fastBandPass(LFP(a:z,i),1e4,10);
	catch
	end
end

figure('outerposition',[0 0 900 700],'PaperUnits','points','PaperSize',[900 700]); hold on
clear a
for i = 1:3
	a(i) = subplot(3,1,i); hold on
	set(a(i),'XLim',[20 50])
end
time = 1e-3*(1:length(PID));
c = parula(length(unique(orn))+1);
for i = 1:length(unique(orn))
	example_orn = i;


	S = PID(:,orn==example_orn); 
	X = filteredLFP(:,orn==example_orn);
	R = fA(:,orn==example_orn);
	R(:,sum(R) == 0) = [];

	axes(a(1))
	errorShade(time,nanmean(S,2),sem(S'),'Color',c(i,:),'SubSample',50);
	axes(a(2))
	errorShade(time,nanmean(X,2),sem(X'),'Color',c(i,:),'SubSample',50);
	axes(a(3))
	errorShade(time,nanmean(fA(:,orn==example_orn),2),sem(fA(:,orn==example_orn)'),'Color',c(i,:),'SubSample',50);
	
end
ylabel(a(1),'Stimulus (V)')
ylabel(a(2),'\DeltaLFP (mV)')
ylabel(a(3),'Firing Rate (Hz)')
 set(a(2),'YLim',[-4 2.5])

prettyFig;

if being_published
	snapnow
	delete(gcf)
end

%% Summary of all Filters
% We now extract all possible filters from this dataset, for each neuron: 

% K -- PID -> LFP filter
if ~exist('K','var')
	K = NaN(1e3,length(unique(orn)));
	for i = 1:length(unique(orn))
		resp = nanmean(filteredLFP(20e3:55e3,orn==i),2);
		stim = nanmean(PID(20e3:55e3,orn==i),2);
		temp = fitFilter2Data(stim,resp,'reg',1,'filter_length',1400,'offset',400);
		% throw out 200ms on either end
		temp(1:200) = [];
		temp(end-199:end) = [];
		K(:,i) = temp;
	end
end

if ~exist('K2','var')
	K2 = NaN(1e3,length(unique(orn)));
	for i = 1:length(unique(orn))
		stim = nanmean(filteredLFP(20e3:55e3,orn==i),2);
		resp = nanmean(fA(20e3:55e3,orn==i),2);
		temp = fitFilter2Data(stim,resp,'reg',1,'filter_length',1400,'offset',400);
		% throw out 200ms on either end
		temp(1:200) = [];
		temp(end-199:end) = [];
		K2(:,i) = temp;
	end
end

if ~exist('K3','var')
	K3 = NaN(1e3,length(unique(orn)));
	for i = 1:length(unique(orn))
		stim = nanmean(PID(20e3:55e3,orn==i),2);
		resp = nanmean(fA(20e3:55e3,orn==i),2);
		temp = fitFilter2Data(stim,resp,'reg',1,'filter_length',1400,'offset',400);
		% throw out 200ms on either end
		temp(1:200) = [];
		temp(end-199:end) = [];
		K3(:,i) = temp;
	end
end

figure('outerposition',[0 0 1400 450],'PaperUnits','points','PaperSize',[1400 550]); hold on
clear a
for i = 1:3
	a(i) = subplot(1,3,i); hold on
	xlabel('Lag (s)')
end
filtertime = 1e-3*(1:length(K)) - .2;
for i = 1:width(K)
	axes(a(1))
	plot(filtertime,K)
	axes(a(2))
	plot(filtertime,K2)
	axes(a(3))
	plot(filtertime,K3)

end
title(a(1),'PID \rightarrow LFP Filter')
title(a(2),'LFP \rightarrow Firing rate Filter')
title(a(3),'PID \rightarrow Firing rate Filter')

prettyFig;

if being_published
	snapnow
	delete(gcf)
end

%% Whiff Analysis
% I now analyse the data like I analysed the naturalistic stimulus data -- by whiff. I also fit a Hill function to the LFP responses to see what the $n$ parameter is. 

clear data

for i = 1:max(orn);
	data(i).S = nanmean(PID(20e3:50e3,orn==i),2); 
	data(i).X = nanmean(filteredLFP(20e3:50e3,orn==i),2);
	R = fA(20e3:50e3,orn==i);
	R(:,sum(R) == 0) = [];
	data(i).R = nanmean(R,2);


end

ft = fittype(' hillFit(x,A,k,n,x_offset)');

figure('outerposition',[0 0 1400 611],'PaperUnits','points','PaperSize',[1400 611]); hold on
for i = 1:length(data)
	S = data(i).S;
	X = data(i).X;
	R = data(i).R;
	ws = whiffStatistics(S,X,R,300,'MinPeakProminence',max(S/1e2),'debug',false);

	subplot(2,length(data),i); hold on
	x = ws.stim_peaks;
	y = -ws.peak_LFP; y = y  - min(y);
	plot(x,y,'.','MarkerSize',20,'Color','k')
	set(gca,'XScale','log')
	ylabel('LFP_{peak} (mV)')
	axis square
	% fit a Hill function to this
	ff = fit(x(~isnan(y)),y(~isnan(y)),ft,'StartPoint',[25 .1 2 0],'Lower',[1 1e-3 .1 -10]);
	l = plot(sort(x),ff(sort(x)),'r');
	legend(l,['n = ' oval(ff.n)],'Location','southeast')

	subplot(2,length(data),length(data)+i); hold on
	x = ws.stim_peaks; 
	y = ws.peak_firing_rate;
	plot(x,y,'.','MarkerSize',20,'Color','k')
	% fit a Hill function to this
	ff = fit(x(~isnan(y)),y(~isnan(y)),ft,'StartPoint',[25 .1 2 0],'Lower',[1 1e-3 .1 -10]);
	l = plot(sort(x),ff(sort(x)),'r');
	legend(l,['n = ' oval(ff.n)],'Location','southeast')

	set(gca,'XScale','log','YLim',[0 70])
	xlabel('Stim_{peak} (V)')
	ylabel('Firing rate_{peak} (mV)')
	axis square

end

prettyFig();

if being_published
	snapnow
	delete(gcf)
end


%%
% It looks like $n<1$. However, this estimate complicated by the fact that I don't know what the absolute value of the LFP is (because of the way the experiment was performed). If I introduce an offset in the y-axis, I can get a very different $n$. Is there another way to estimate $n$? 

%% Is there fast gain control? 
% Now I collect whiffs with the same amplitude, and look at the LFP and firing rate responses to see if there is any difference that could be attributed to stimulus context. 


ws = plotScaledNatStimWhiffStats(data,false);
s_range = [.313 .329];

figure('outerposition',[0 0 1501 502],'PaperUnits','points','PaperSize',[1501 502]); hold on
for i = 1:3
	ax(i) = subplot(1,3,i); hold on
end

show_these = ws(3).stim_peak_loc((ws(3).stim_peaks>s_range(1) & ws(3).stim_peaks<s_range(2)));


for i = 3
	for j = 1:length(show_these)
		this_loc = show_these(j);

		S = data(i).S;
		X = data(i).X;
		R = data(i).R;

		a = this_loc - 300;
		z = this_loc+300;

		if a > 1 && z < length(S)
			plot(ax(1),S(a:z))
			plot(ax(2),X(a:z))
			plot(ax(3),R(a:z))
		end

	end
end
ylabel(ax(1),'Stimulus (V)')
ylabel(ax(2),'LFP (mV)')
ylabel(ax(3),'Firing rate (Hz)')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%% Fitting an adapting NL model to the data
% In this section, I fit an adapting NL and NLN model to the LFP and firing rate to determine what the $n$ is in this data, and to see if that model can account for this data. 

% clear fd
% for i = 1:max(orn)
% 	fd(i).stimulus = data(i).S;
% 	fd(i).response = data(i).R;
% 	fd(i).response(1:1e3) = NaN;
% end

% for i = 1:max(orn)
% 	p(i) = fitModel2Data(@aNLN5,fd(i),'nsteps',50,'use_parallel',true,'p0',p(i));
% end


%%
% In the following figure, I show the predictions made by adapting NL and NLN models for the LFP and the firing rate. Models were fit individually to each neuron. Each column shows the model fits and data from a single neuron. 

load('aNL5_fits_to_LVF.mat','p')
p_LFP = p;
load('aNLN5_fits_to_LVF.mat','p')


figure('outerposition',[0 0 1550 811],'PaperUnits','points','PaperSize',[1550 811]); hold on
for i = 1:max(orn)

	% show firing rate fits
	subplot(2,max(orn),i); hold on
	R = aNL5(data(i).S,p_LFP(i));
	plot(-R(1e3:end),data(i).X(1e3:end),'k.')
	legend(['r^2 = ' oval(rsquare(-R(1e3:end),data(i).X(1e3:end)))],'Location','southeast')
	xlabel('NL model prediction (mV)')
	ylabel('\Delta LFP (mV)')

	% show firing rate fits
	subplot(2,max(orn),max(orn)+i); hold on
	R = aNLN5(data(i).S,p(i));
	plot(R(1e3:end),data(i).R(1e3:end),'k.')
	legend(['r^2 = ' oval(rsquare(R(1e3:end),data(i).R(1e3:end)))],'Location','southeast')
	xlabel('NLN model prediction (Hz)')
	ylabel('Firing rate (Hz)')
end

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% What are the $n$ parameters in these fits?

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(1:max(orn),[p.Hill_n],'k+-')
plot(1:max(orn),[p_LFP.Hill_n],'r+-')
legend({'Firing rate','LFP'})
xlabel('ORN #')
ylabel('n in fit')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%% Timescale of gain control. 
% In this section, I attempt to estimate the timescale of gain control from model fits to this data. Since in this data, there are large, variable excursions on a background, we expect gain control to play an important role. In this model, $k_D$ varies accoridng to:
% 
% $$\dot{k}_{D}=B*k_{D}\left(a-\nicefrac{1}{2}\right)$$
% 

%%
% where B is a free parameter that controls both the "amount" of gain control and its timescale. 



%%
% In the following figure, in each plot, each line/colour corresponds to data from a single neuron. (a) Shows the stimulus distribution. Note that it is on a log scale. (b) Shows the distribution of $k_D$ around its mean. In the model fits, there appears to be some gain control that changes $k_D$ over the course of the stimulus. How significant is this? (c) shows the ratio of standard deviations of $k_D$ to the standard deviation in the stimulus. Every neuron exhibits non-zero gain control, in the framework of an adapting NLN model. Now, how quickly does the $k_D# vary? (d) shows the autocorrelation functions of $k_D$ for each neuron, cf. autocorrelation functions for the corresponding stimulus. In (e-f), I vary the parameter $B$ in the model, which governs the timescale of gain control (and how much of it is there), and plot $r^2$ of the model vs. the data. In (f), I plot the $r^2$ of the model prediction vs. the autocorrelation time of $k_D$. 


%%
% Now I make the same plots, but for the firing rate models. 

figure('outerposition',[0 0 1333 801],'PaperUnits','points','PaperSize',[1333 801]); hold on
clear ax
for i = 1:6
	ax(i) = subplot(2,3,i); hold on
end

c = lines(max(orn));
for i = 1:max(orn)

	S = nanmean(PID(:,orn==i),2); 
	R = nanmean(fA(:,orn==i),2);
	R(:,sum(R)==0)= [];
	R = nanmean(R,2);


	[hy,hx] = histcounts(S);
	hy = hy/sum(hy);
	hx = hx(2:end) + mean(diff(hx));
	plot(ax(1),hx,hy,'Color',c(i,:))

	% run it through the model
	[~, ~, ~, kD] = aNLN5(S,p(i));
	plot(ax(3),i,std(kD)/std(S),'+','Color',c(i,:),'MarkerSize',20)

	[hy,hx] = histcounts(kD/mean(kD));
	hy = hy/sum(hy);
	hx = hx(2:end) + mean(diff(hx));
	plot(ax(2),hx,hy,'Color',c(i,:))

	% compute autocorrelation of KD
	kD = kD(10e3:end);
	kD = kD - mean(kD);
	kD = kD/std(kD);
	[a,lags] = autocorr(kD,length(kD)-1);
	plot(ax(4),lags,a,':','Color',c(i,:),'LineWidth',2)

	stim = S(10e3:end);
	stim = stim - mean(stim);
	stim = stim/std(stim);
	[a,lags] = autocorr(stim,length(stim)-1);
	plot(ax(4),lags,a,'Color',c(i,:),'LineWidth',2)

	% now, vary B and calcualte r^2 and autocorrelation of KD
	all_B = logspace(log10(p(i).B/10),log10(p(i).B*10),21);
	r2 = NaN*all_B;
	tau = NaN*all_B;
	q = p(i);
	for j = 1:length(all_B)
		q.B = all_B(j);
		[Rp, ~, ~, kD] = aNLN5(S,q);
		r2(j) = rsquare(Rp(20e3:50e3),R(20e3:50e3));
		% also compute the autocorrelation time
		kD = kD(10e3:end);
		kD = kD - mean(kD);
		kD = kD/std(kD);
		[a,lags] = autocorr(kD,length(kD)-1);
		tau(j) = lags(find(a<1/exp(1),1,'first'));
	end
	x = logspace(-1,1,21);
	plot(ax(5),x,r2,'Color',c(i,:))
	plot(ax(6),tau,r2,'Color',c(i,:))
end

xlabel(ax(1),'Stimulus (V)')
ylabel(ax(1),'Probability')
set(ax(1),'XScale','log')

xlabel(ax(2),'k_{D}/<k_D>')
ylabel(ax(2),'Probability')
set(ax(2),'XScale','linear','XLim',[.1 2])

set(ax(3),'YLim',[0 .4])
xlabel(ax(3),'ORN #')
ylabel(ax(3),'\sigma_{k_D}/\sigma_{S}')

clear l L
l(1) = plot(ax(4),NaN,NaN,'k:');
l(2) = plot(ax(4),NaN,NaN,'k','LineWidth',2);
legend(l,{'k_D','Stimulus'},'Location','southwest')

set(ax(4),'XScale','log','XLim',[10 1e4],'XTick',[1 10 100 1e3 1e4])
xlabel(ax(4),'Lag (ms)')
ylabel(ax(4),'autocorrelation')

xlabel(ax(5),'B/B_{fit}')
set(ax(5),'XScale','log','XLim',[.1 10])
ylabel(ax(5),'r^2')

xlabel(ax(6),'\tau_{k_D} (ms)')
set(ax(6),'XScale','log','XLim',[100 1e4])
ylabel(ax(6),'r^2')

prettyFig();

labelFigure

if being_published
	snapnow
	delete(gcf)
end

%%
% This looks promising. The fundamental problem here is that in this model, the amount of gain control and the timescale of gain control are deeply linked. This is a problem because we don't know if the model fits the data well because it's getting the timescale of gain control right, or because it's getting the amount right. This method of analysis only tells us what the best value of the parameter B that reproduces the data is -- going from that to the timescale is a bit of a jump.

%%
% What we instead need is a model where the amount of gain control and the timescale of gain control are seperated. 


% clear fd
% for i = 1:max(orn)
% 	fd(i).stimulus = data(i).S;
% 	fd(i).response = data(i).R;
% 	fd(i).response(1:1e3) = NaN;
% end

% for i = max(orn):-1:1
% 	p(i) = fitModel2Data(@aNLN3,fd(i),'nsteps',50,'use_parallel',true,'p0',p(i));
% end

%%
% In the following figure, I plot the data vs. the predictions of this model to see how good a fit we can get. 



load('aNLN3_fits_to_LVF.mat','p')



figure('outerposition',[0 0 1550 501],'PaperUnits','points','PaperSize',[1550 501]); hold on
for i = 1:max(orn)
	% show firing rate fits
	subplot(1,max(orn),i); hold on
	R = aNLN3(data(i).S,p(i));
	plot(R(1e3:end),data(i).R(1e3:end),'k.')
	legend(['r^2 = ' oval(rsquare(R(1e3:end),data(i).R(1e3:end)))],'Location','southeast')
	xlabel('NLN model prediction (Hz)')
	ylabel('Firing rate (Hz)')
end

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% It seems really good. However, does the model actually use the front-end gain control to achieve this impressive fit? To determine this, I plot the time-dependent $k_D$ as a function of time for the various neuron fits. 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
for i = 1:max(orn)

	S = nanmean(PID(:,orn==i),2); 
	[~, ~, ~, kD] = aNLN3(S,p(i));
	plot(time,kD)

end
xlabel('Time (s)')
ylabel('k_D')
set(gca,'XLim',[0 60],'YLim',[0 0.3])

prettyFig();

if being_published
	snapnow
	delete(gcf)
end


%%
% It looks like the gain control mechanism isn't being involved at all. Since the fits are really good, I conclude that there is no gain control going on in this data. 



%% Version Info

pFooter;