% Fold Change Analysis
% 

pHeader;

%% Fold Change Analysis
% In this document, we analyse data from an experiment to determine if ORNs do fold change detection. 

%% The Stimulus
% To study fold change detection, we first need to deliver a stimulus that has a fold change. Basically, we want a step of stimulus, such that the post-step stimulus is twice the pre-step stimulus, for various background levels. However, this proved to be very hard to engineer. Instead, we went for a shotgun approach, where steps of various heights were presented on top of a background. By chance, we will get a bunch of them to have a fold change of close to 2. The following figure shows the result of this approach. 

% get the data

[PID, LFP, fA] = consolidateData('/local-data/DA-paper/fold-change/ab3',1);

rm_this = sum(fA) == 0;
fA(:,rm_this) = [];
PID(:,rm_this) = [];
LFP(:,rm_this) = [];


% remove the baseline from the LFP
for i = 1:width(LFP)
	LFP(:,i) = LFP(:,i) - mean(LFP(1:5e3,i));
end

% differentiate the LFP
dLFP = LFP;
for i = 1:width(PID)
	% first high pass them to remove spikes
	LFP(:,i) = bandPass(LFP(:,i),Inf,30);
	dLFP(:,i) = -1e4*filtfilt(ones(10,1),10,[0; diff(LFP(:,i))]);
end

fold_change = mean(PID(10.5e3:11e3,:))./mean(PID(6e3:10e3,:));
background_stim = mean(PID(6e3:10e3,:));

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(fold_change,background_stim,'k+')
plot([1.8 1.8],[.1 .6],'k--')
plot([2.2 2.2],[.1 .6],'k--')
xlabel('Fold Change')
ylabel('Background Stimulus')

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end


%%
% We now plot the stimulus normalised to the pre-step level to show that our fold-change stimulus is reasonable. In the following figure, we plot each trial seperately, and colour it by the pre-step stimulus. 

rm_this = fold_change > 2.2 | fold_change < 1.8;
PID(:,rm_this) = [];
fA(:,rm_this) = [];
LFP(:,rm_this) = [];
dLFP(:,rm_this) = [];

fold_change = mean(PID(10.5e3:11e3,:))./mean(PID(6e3:10e3,:));
background_stim = mean(PID(6e3:10e3,:));


c = parula(100);
cc = background_stim - min(background_stim);
cc = cc/max(cc); cc = 1+round(cc*99);

time = 1e-3*(1:length(PID));
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
for i = 1:width(PID)
	plot(time,PID(:,i),'Color',c(cc(i),:))
end
set(gca,'XLim',[8 14])
xlabel('Time (s)')
ylabel('Stimulus (V)')

subplot(1,2,2), hold on
for i = 1:width(PID)
	x = PID(:,i)/mean(PID(6e3:10e3,i));
	plot(time,x,'Color',c(cc(i),:))
end
set(gca,'XLim',[8 14])
xlabel('Time (s)')
ylabel('Stimulus (norm)')

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%% ORN Response
% We now plot the neuron's response to these fold change stimuli. In the following traces, each trial is plotted separately, and is colour coded by the pre-step stimulus level.


figure('outerposition',[0 0 1000 800],'PaperUnits','points','PaperSize',[1000 800]); hold on 
subplot(2,2,1), hold on
for i = 1:width(PID)
	plot(time,LFP(:,i),'Color',c(cc(i),:))
end
xlabel('Time (s)')
set(gca,'XLim',[8 14])
ylabel('LFP (mV)')

subplot(2,2,2), hold on
for i = 1:width(PID)
	x = LFP(:,i) - mean(LFP(9e3:10e3,i));
	plot(time,x,'Color',c(cc(i),:))
end
xlabel('Time (s)')
set(gca,'XLim',[8 14])
ylabel('\DeltaLFP (mV)')

subplot(2,2,3), hold on
for i = 1:width(PID)
	plot(time,dLFP(:,i),'Color',c(cc(i),:))
end
xlabel('Time (s)')
set(gca,'XLim',[8 14])
ylabel('dLFP (mV/s)')

subplot(2,2,4), hold on
for i = 1:width(PID)
	plot(time,fA(:,i),'Color',c(cc(i),:))
end
xlabel('Time (s)')
set(gca,'XLim',[8 14])
ylabel('Firing Rate (Hz)')

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%% Version Info
%
pFooter;

