% obp28a_analysis_ab8.m
% 
% created by Srinivas Gorur-Shandilya at 3:42 , 01 January 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.ts


%% OBP28a analysis of Gain
% In this document, we see if the gain properties of the ab8A ORN are changed when we delete the OBP 28a. 


pHeader;

% prep and merge data
use_cache = 1;
[PID, LFP, fA, paradigm, orn,fly] = consolidateData('/local-data/obp/ab8/ko',use_cache);
geno = ones(length(orn),1);

[PID2, LFP2, fA2, paradigm2, orn2,fly2] = consolidateData('/local-data/obp/ab8/wcs',use_cache);
geno2 = 2*ones(length(orn2),1);

PID = [PID PID2]; clear PID2
fA = [fA fA2]; clear fA2
LFP = [LFP LFP2]; clear LFP2
paradigm = [paradigm paradigm2]; clear paradigm2
orn = [orn orn2]; clear orn2
fly = [fly; fly2]; clear fly2
geno = [geno; geno2]; clear geno2
paradigm = paradigm(:);
orn = orn(:);


% clean up data
rm_this = isnan(sum(LFP));
fA(:,rm_this) = [];
PID(:,rm_this) = [];
LFP(:,rm_this) = [];
paradigm(rm_this) = [];
orn(rm_this) = [];
fly(rm_this) = [];
geno(rm_this) = [];
fA(:,sum(fA) == 0) = NaN;


% bandPass LFP
filtered_LFP = LFP;
for i = 1:width(LFP)
	filtered_LFP(:,i) = bandPass(LFP(:,i),1000,10);
end

a = 10e3; z = 50e3;
[K1,LFP_pred,LFP_gain] = extractFilters(PID,filtered_LFP,'use_cache',true,'a',a,'z',z);

[K2,fA_pred,fA_gain] = extractFilters(PID,fA,'use_cache',true,'a',a,'z',z);

%% Overview of the data
% In the following figure, we plot all LFP and firing rate traces for visual inspection. 

figure('outerposition',[0 0 1000 800],'PaperUnits','points','PaperSize',[1000 800]); hold on
c = parula(max(paradigm));
subplot(3,1,1), hold on
time = 1e-3*(1:length(PID));
for i = 1:max(paradigm)-1
	plot(time,nanmean(PID(:,paradigm == i),2),'Color',c(i,:));
end
ylabel('PID (V)')
set(gca,'XLim',[0 60])

subplot(3,1,2), hold on
for i = 1:max(paradigm)-1
	plot(time,nanmean(filtered_LFP(:,paradigm == i),2),'Color',c(i,:));
end
ylabel('LFP (mV)')
set(gca,'XLim',[0 60])

subplot(3,1,3), hold on
for i = 1:max(paradigm)-1
	plot(time,nanmean(fA(:,paradigm == i),2),'Color',c(i,:));
end
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')
set(gca,'XLim',[0 60])

prettyFig('fs=20;')

if being_published
	snapnow
	delete(gcf)
end

%%
% In the following figure, we compare the KO to the control genotype.

figure('outerposition',[0 0 1000 800],'PaperUnits','points','PaperSize',[1000 800]); hold on
c = parula(max(paradigm));
subplot(3,1,1), hold on
time = 1e-3*(1:length(PID));
for i = 1:max(paradigm)-1
	plot(time,nanmean(PID(:,paradigm == i & geno == 2),2),'Color','k');
	plot(time,nanmean(PID(:,paradigm == i & geno == 1),2),'Color','r');
end
ylabel('PID (V)')
set(gca,'XLim',[30 40])

subplot(3,1,2), hold on
for i = 1:max(paradigm)-1
	plot(time,nanmean(filtered_LFP(:,paradigm == i & geno == 2),2),'Color','k');
	plot(time,nanmean(filtered_LFP(:,paradigm == i & geno == 1),2),'Color','r');
end
ylabel('LFP (mV)')
set(gca,'XLim',[30 40])

subplot(3,1,3), hold on
for i = 1:max(paradigm)-1
	plot(time,nanmean(fA(:,paradigm == i & geno == 2),2),'Color','k');
	plot(time,nanmean(fA(:,paradigm == i & geno == 1),2),'Color','r');
end
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')
set(gca,'XLim',[30 40])

prettyFig('fs=20;')

if being_published
	snapnow
	delete(gcf)
end

%% Filters and Output Nonlinearities 
% In this section, we compare the filter shapes and the output nonlinearities for each dose, for each odour dose concentration. The figure below shows this comparison for the PID to LFP transformation. 

figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[900 800]); hold on
filtertime = 1e-3*(1:length(K1)) - .1;
opacity = .3;
for i = 1:max(paradigm)-1
	subplot(3,2,2*(i-1)+1), hold on
	this_K = K1(:,paradigm==i & geno == 2);
	errorShade(filtertime,nanmean(this_K,2),nanstd(this_K'),'Color','k');
	this_K = K1(:,paradigm==i & geno == 1);
	errorShade(filtertime,nanmean(this_K,2),nanstd(this_K'),'Color','r');
	xlabel('Filter Lag (s)')

	subplot(3,2,2*(i-1)+2), hold on
	x = LFP_pred(a:z,paradigm==i & geno == 2);
	y = filtered_LFP(a:z,paradigm==i & geno == 2);
	for j = 1:width(x)
		try
			[~,data]=plotPieceWiseLinear(x(:,j),y(:,j),'nbins',50,'make_plot',false);
			h = plot(data.x,data.y);
			h.Color = [0 0 0 opacity];
		end
	end
	x = LFP_pred(a:z,paradigm==i & geno == 1);
	y = filtered_LFP(a:z,paradigm==i & geno == 1);
	for j = 1:width(x)
		try
			[~,data]=plotPieceWiseLinear(x(:,j),y(:,j),'nbins',50,'make_plot',false);
			h = plot(data.x,data.y);
			h.Color = [1 0 0 opacity];
		end
	end
	xlabel('Linear Prediction (mV)')
	ylabel('\DeltaLFP(mV)')
end

prettyFig('fs=20;')

if being_published
	snapnow
	delete(gcf)
end

%%
% Now we compare filters and nonlinearities for the PID->firing rate transformation. 

figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[900 800]); hold on
filtertime = 1e-3*(1:length(K2)) - .1;
for i = 1:max(paradigm)-1
	subplot(3,2,2*(i-1)+1), hold on
	this_K = K2(:,paradigm==i & geno == 2);
	errorShade(filtertime,nanmean(this_K,2),nanstd(this_K'),'Color','k');
	this_K = K2(:,paradigm==i & geno == 1);
	errorShade(filtertime,nanmean(this_K,2),nanstd(this_K'),'Color','r');
	xlabel('Filter Lag (s)')

	subplot(3,2,2*(i-1)+2), hold on
	x = fA_pred(a:z,paradigm==i & geno == 2);
	y = fA(a:z,paradigm==i & geno == 2);
	for j = 1:width(x)
		try
			[~,data]=plotPieceWiseLinear(x(:,j),y(:,j),'nbins',50,'make_plot',false);
			h = plot(data.x,data.y);
			h.Color = [0 0 0 opacity];
		end
	end
	x = fA_pred(a:z,paradigm==i & geno == 1);
	y = fA(a:z,paradigm==i & geno == 1);
	for j = 1:width(x)
		try
			[~,data]=plotPieceWiseLinear(x(:,j),y(:,j),'nbins',50,'make_plot',false);
			h = plot(data.x,data.y);
			h.Color = [1 0 0 opacity];
		end
	end
	xlabel('Linear Prediction (V)')
	ylabel('Firing Rate(Hz)')
end

prettyFig('fs=20;')

if being_published
	snapnow
	delete(gcf)
end


%% Comparison of Gain 
% In the following figure we compare the gain in the LFP and in the firing rate between the wCS control (black) and the CRISPR OBP deletion (red). 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
x = mean(PID(a:z,geno==1 & paradigm < 4)); y = LFP_gain(geno==1 & paradigm < 4);
ff = fit(x(:),y(:),'power1');
plot(.1:0.1:2,ff(.1:0.1:2),'r')
plot(x,y,'r+');

x = mean(PID(a:z,geno==2 & paradigm < 4)); y = LFP_gain(geno==2 & paradigm < 4 );
ff = fit(x(:),y(:),'power1');
plot(.1:0.1:2,ff(.1:0.1:2),'k')
plot(x,y,'k+');

set(gca,'XScale','log','YScale','log','XLim',[.1 10])
xlabel('Mean Stimulus (V)')
ylabel('LFP Gain (mV/V)')

subplot(1,2,2), hold on
x = mean(PID(a:z,geno==1 & paradigm < 4)); y = fA_gain(geno==1 & paradigm < 4);
x(isnan(y)) = []; y(isnan(y)) = [];
ff = fit(x(:),y(:),'power1');
plot(.1:0.1:2,ff(.1:0.1:2),'r')
plot(x,y,'r+');

x = mean(PID(a:z,geno==2 & paradigm < 4)); y = fA_gain(geno==2 & paradigm < 4);
x(isnan(y)) = []; y(isnan(y)) = [];
ff = fit(x(:),y(:),'power1');
plot(.1:0.1:2,ff(.1:0.1:2),'k')
plot(x,y,'k+');

set(gca,'XScale','log','YScale','log','XLim',[.1 10])
xlabel('Mean Stimulus (V)')
ylabel('Firing Gain (Hz/V)')


prettyFig('fs=20;')

if being_published
	snapnow
	delete(gcf)
end

%% Ablation of the B neuron
% In an effort to make the spikes more sortable, Nikki ablated the B neuron usind DsArm. Unfortunately, for unknown reasons, the B neuron was not ablated. In the following figure, we try to compare responses from the OBP KO and the wCS control, where we hoped the B neuron was ablated, but in reality, it wasn't. 


% prep and merge data
use_cache = 1;
[PID, LFP, fA, paradigm, orn,fly] = consolidateData('/local-data/obp/ab8/B-ablated-ko',use_cache);
geno = ones(length(orn),1);

[PID2, LFP2, fA2, paradigm2, orn2,fly2] = consolidateData('/local-data/obp/ab8/B-ablated-wcs',use_cache);
geno2 = 2*ones(length(orn2),1);

PID = [PID PID2]; clear PID2
fA = [fA fA2]; clear fA2
LFP = [LFP LFP2]; clear LFP2
paradigm = [paradigm paradigm2]; clear paradigm2
orn = [orn orn2]; clear orn2
fly = [fly; fly2]; clear fly2
geno = [geno; geno2]; clear geno2
paradigm = paradigm(:);
orn = orn(:);

% clean up data
rm_this = isnan(sum(LFP));
fA(:,rm_this) = [];
PID(:,rm_this) = [];
LFP(:,rm_this) = [];
paradigm(rm_this) = [];
orn(rm_this) = [];
fly(rm_this) = [];
geno(rm_this) = [];
fA(:,sum(fA) == 0) = NaN;


% bandPass LFP
filtered_LFP = LFP;
for i = 1:width(LFP)
	filtered_LFP(:,i) = bandPass(LFP(:,i),1000,10);
end

a = 10e3; z = 50e3;
[K1,LFP_pred,LFP_gain] = extractFilters(PID,filtered_LFP,'use_cache',true,'a',a,'z',z);

[K2,fA_pred,fA_gain] = extractFilters(PID,fA,'use_cache',true,'a',a,'z',z);

%%
% In the following figure, we compare filters and nonlinearities for the stimulus to LFP transformations. 

figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[900 800]); hold on
filtertime = 1e-3*(1:length(K1)) - .1;
opacity = .3;
for i = 1:max(paradigm)
	subplot(3,2,2*(i-1)+1), hold on
	this_K = K1(:,paradigm==i & geno == 2);
	errorShade(filtertime,nanmean(this_K,2),nanstd(this_K'),'Color','k');
	this_K = K1(:,paradigm==i & geno == 1);
	errorShade(filtertime,nanmean(this_K,2),nanstd(this_K'),'Color','r');
	xlabel('Filter Lag (s)')

	subplot(3,2,2*(i-1)+2), hold on
	x = LFP_pred(a:z,paradigm==i & geno == 2);
	y = filtered_LFP(a:z,paradigm==i & geno == 2);
	for j = 1:width(x)
		try
			[~,data]=plotPieceWiseLinear(x(:,j),y(:,j),'nbins',50,'make_plot',false);
			h = plot(data.x,data.y);
			h.Color = [0 0 0 opacity];
		end
	end
	x = LFP_pred(a:z,paradigm==i & geno == 1);
	y = filtered_LFP(a:z,paradigm==i & geno == 1);
	for j = 1:width(x)
		try
			[~,data]=plotPieceWiseLinear(x(:,j),y(:,j),'nbins',50,'make_plot',false);
			h = plot(data.x,data.y);
			h.Color = [1 0 0 opacity];
		end
	end
	xlabel('Linear Prediction (mV)')
	ylabel('\DeltaLFP(mV)')
end

prettyFig('fs=20;')

if being_published
	snapnow
	delete(gcf)
end

%%
% Now we compare filters and nonlinearities for the PID->firing rate transformation. 

figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[900 800]); hold on
filtertime = 1e-3*(1:length(K2)) - .1;
for i = 1:max(paradigm)
	subplot(3,2,2*(i-1)+1), hold on
	this_K = K2(:,paradigm==i & geno == 2);
	errorShade(filtertime,nanmean(this_K,2),nanstd(this_K'),'Color','k');
	this_K = K2(:,paradigm==i & geno == 1);
	errorShade(filtertime,nanmean(this_K,2),nanstd(this_K'),'Color','r');
	xlabel('Filter Lag (s)')

	subplot(3,2,2*(i-1)+2), hold on
	x = fA_pred(a:z,paradigm==i & geno == 2);
	y = fA(a:z,paradigm==i & geno == 2);
	for j = 1:width(x)
		try
			[~,data]=plotPieceWiseLinear(x(:,j),y(:,j),'nbins',50,'make_plot',false);
			h = plot(data.x,data.y);
			h.Color = [0 0 0 opacity];
		end
	end
	x = fA_pred(a:z,paradigm==i & geno == 1);
	y = fA(a:z,paradigm==i & geno == 1);
	for j = 1:width(x)
		try
			[~,data]=plotPieceWiseLinear(x(:,j),y(:,j),'nbins',50,'make_plot',false);
			h = plot(data.x,data.y);
			h.Color = [1 0 0 opacity];
		end
	end
	xlabel('Linear Prediction (V)')
	ylabel('Firing Rate(Hz)')
end

prettyFig('fs=20;')


%% Version Info
%

pFooter;