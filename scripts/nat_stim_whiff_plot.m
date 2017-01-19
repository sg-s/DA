pHeader;

%% Detailed look at naturalistic stimulus whiffs
% In this document, I look at the ORN responses (both LFP and firing rate) to individual whiffs in the naturalistic stimulus response, to try to understand what's going on in the response without any filter analysis. 



% get the naturalistic stimuli 
load('/local-data/DA-paper/data-for-paper/nat-stim/ab3A_nat_stim.ORNData','-mat')

% remove baseline from stimulus, LFP
for i = 1:length(od)
	for j = 1:od(i).n_trials
		od(i).stimulus(:,j) = od(i).stimulus(:,j) - min(od(i).stimulus(1:5e3,j));
		od(i).LFP(:,j) = od(i).LFP(:,j) - mean(od(i).LFP(1:5e3,j));
	end
end

% pick an example orn
example_orn = 3;
R = nanmean(od(example_orn).firing_rate,2);
S = nanmean(od(example_orn).stimulus,2);
X = nanmean(od(example_orn).LFP,2);
tA = 1e-3*(1:length(S));

figure('outerposition',[0 0 1601 855],'PaperUnits','points','PaperSize',[1601 855]); hold on

a = [36.36 47.9; 47.9 59.46];
before = .5;
after = .5;

c = lines(2);

% show the stimulus
subplot(3,6,1:3); hold on
plot(tA,S,'k')
ylabel('Stimulus (V)')
set(gca,'XLim',[0 70],'YLim',[-.2 9])

% indicate the segments we zoom into we use
plot([a(1,1)-before a(1,1)+after],[2 2],'Color',c(1,:),'LineWidth',4)
plot([a(1,2)-before a(1,2)+after],[2 2],'Color',c(2,:),'LineWidth',4)

plot([a(2,1)-before a(2,1)+after],[1.5 1.5],'Color',c(1,:),'LineWidth',4)
plot([a(2,2)-before a(2,2)+after],[1.5 1.5],'Color',c(2,:),'LineWidth',4)

% show LFP response
subplot(3,6,7:9); hold on
plot(tA,X,'k')
ylabel('ab3 LFP (\DeltamV)')
set(gca,'XLim',[0 70],'YLim',[-20 5])

% % indicate the segments we zoom into we use
% plot([a(1,1)-before a(1,1)+after],[125 125],'Color',c(1,:),'LineWidth',4)
% plot([a(1,2)-before a(1,2)+after],[125 125],'Color',c(2,:),'LineWidth',4)

% plot([a(2,1)-before a(2,1)+after],[125 125],'Color',c(1,:),'LineWidth',4)
% plot([a(2,2)-before a(2,2)+after],[125 125],'Color',c(2,:),'LineWidth',4)

% show the firing rate
subplot(3,6,13:15); hold on
plot(tA,R,'k')
ylabel('ab3A firing rate (Hz)')
set(gca,'XLim',[0 70],'YLim',[-5 170])

% % indicate the segments we zoom into we use
% plot([a(1,1)-before a(1,1)+after],[125 125],'Color',c(1,:),'LineWidth',4)
% plot([a(1,2)-before a(1,2)+after],[125 125],'Color',c(2,:),'LineWidth',4)

% plot([a(2,1)-before a(2,1)+after],[125 125],'Color',c(1,:),'LineWidth',4)
% plot([a(2,2)-before a(2,2)+after],[125 125],'Color',c(2,:),'LineWidth',4)


stim_plots = [4 5];
lfp_plots = [10 11];
firing_plots = [16 17];
for i = 1:2
	for j = 1:2
		% plot the stimulus
		subplot(3,6,stim_plots(i)); hold on
		s = 1e3*(a(i,j) -  before);
		e = 1e3*(a(i,j) + after);
		x = tA(s:e); x = x - a(i,j);
		y = S(s:e);

		if j == 1
			plot(x,y,'Color',c(i,:))
		else
			plot(x,y,':','Color',c(i,:))
		end
		set(gca,'XLim',[-.5 .5])

		% plot the LFP
		subplot(3,6,lfp_plots(i)); hold on
		s = 1e3*(a(i,j) -  before);
		e = 1e3*(a(i,j) + after);
		x = tA(s:e); x = x - a(i,j);
		y = X(s:e); 

		if i == 1
			y = y - y(1e3*before);
		end

		if j == 1
			plot(x,y,'Color',c(i,:))
		else
			plot(x,y,':','Color',c(i,:))
		end
		set(gca,'XLim',[-.5 .5])

		% plot the firing rate
		subplot(3,6,firing_plots(i)); hold on
		s = 1e3*(a(i,j) -  before);
		e = 1e3*(a(i,j) + after);
		x = tA(s:e); x = x - a(i,j);
		y = R(s:e);

		if j == 1
			plot(x,y,'Color',c(i,:))
		else
			plot(x,y,':','Color',c(i,:))
		end
		set(gca,'XLim',[-.5 .5])
	end
end

subplot(3,6,4); hold on
set(gca,'YLim',[-.1 .7])

subplot(3,6,5); hold on
set(gca,'YLim',[-.1 .7])

subplot(3,6,13:15); hold on
xlabel('Time (s)')


prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%% Are responses saturated?
% The big question in this dataset is if the responses here are always saturated, and if the stimulus we used is simply too high. We have in this dataset one experiment where the stimulus we used was much smaller. Does this mean that the responses are smaller too? 

%%
% In the following figure, I compare the stimulus and response for two neurons (ab3A and ab2A), one where we used a very small stimulus, and another where we used the normal stimulus. 

time = 1e-3*(1:length(od(1).stimulus));

figure('outerposition',[0 0 1000 901],'PaperUnits','points','PaperSize',[1000 901]); hold on
subplot(3,1,1); hold on
plot(time,mean(od(2).stimulus,2),'r')
plot(time,mean(od(1).stimulus,2),'k')
set(gca,'XLim',[1 68])
ylabel('Stimulus (V)')

subplot(3,1,2); hold on
plot(time,mean(od(2).LFP,2),'r')
plot(time,mean(od(1).LFP,2),'k')
ylabel('\Delta LFP (mV)')
set(gca,'XLim',[1 68])

subplot(3,1,3); hold on
plot(time,mean(od(2).firing_rate,2),'r')
plot(time,mean(od(1).firing_rate,2),'k')
ylabel('Firing rate (Hz)')
xlabel('Time (s)')
set(gca,'XLim',[1 68])

prettyFig();

if being_published
	snapnow
	delete(gcf)
end


%%
% We notice something very interesting. Even though there is a 7-fold difference in the stimulus amplitude, the LFP and firing rate response amplitudes look identical. Even more curiously, the kinetics of the LFP are markedly different: they're much faster in the lower stimulus case. Are the two neurons adapted to different degrees by the two stimuli? To check, I fit a NL model to each dataset. The following figure shows the input nonlinearities of the best-fit models for the two cases. 

for i = 1:2
	data(i).response = mean(od(i).firing_rate,2);
	data(i).stimulus = mean(od(i).stimulus,2);
end

clear p
p(1).   k0 = 0.0157;
p(1).tau_z = 1;
p(1).    B = 0;
p(1).  n_z = 1;
p(1).    n = 1;
p(1). tau1 = 21.6250;
p(1). tau2 = 16.0820;
p(1).  n_y = 4.9492;
p(1).    A = 0.8891;
p(1).    C = 669.2250;


p(2).   k0 = 0.1329;
p(2).tau_z = 1;
p(2).    B = 0;
p(2).  n_z = 1;
p(2).    n = 1;
p(2). tau1 = 11.4688;
p(2). tau2 = 10.6953;
p(2).  n_y = 7.5352;
p(2).    A = 0.8891;
p(2).    C = 892.8285;

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1); hold on
set(gca,'XScale','log','XTick',[1e-3 1e-2 1e-1 1 10])
xlabel('S (V)')
c = [0 0 0; 1 0 0];
for i = 1:2
	S = logspace(-3,1,200);
	x = (S.^p(i).n)./(S.^p(i).n+p(i).k0.^p(i).n);
	plot(S,x,'Color',c(i,:))
end
legend({'Low Stimulus','High Stimulus'},'Location','southeast')

subplot(1,2,2); hold on
xlabel('S/<S>')
set(gca,'XScale','log','XTick',[1e-2 1e-1 1 10 100])
for i = 1:2
	S = logspace(-3,1,200);
	x = (S.^p(i).n)./(S.^p(i).n+p(i).k0.^p(i).n);
	S = S/mean(mean(od(i).stimulus));
	plot(S,x,'Color',c(i,:))
end
plot([1 1],[0 1],'k:')
plot([1e-2 1e2],[0.5 .5],'k:')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% There is something very bizarre going on here. Our crude observation that the neuron seemed to be more sensitive to the smaller stimulus (because its response was identical to the neuron stimulated with the bigger stimulus) is confirmed by the fact that the best-fit input nonlinearities are significantly different. But if we rescale these input nonlinearities by the *mean stimulus presented to the neuron*, the input nonlinearities collapse onto one another. That's weird -- how do these neurons know to do this? Surely it can't be just by chance that these two neurons had the exact $k_D$ that co-incidentally matched the stimulus? 

%% Version Info
%
pFooter;


return

%% Quantification
% In the following section, I quantify some properties of the neuron response. 

% build some metrics for all the data we have
clear whiffs
whiffs.stim_peaks = [];
whiffs.stim_300ms_before_whiff = [];
whiffs.min_stim_before_whiff = [];
whiffs.peak_LFP = [];
whiffs.LFP_before_whiff = [];
whiffs.peak_firing_rate = [];
whiffs.firing_rate_before_whiff = [];
whiffs.orn = [];
for j = 1:length(od)

	R = nanmean(od(j).firing_rate,2);
	S = nanmean(od(j).stimulus,2);
	X = nanmean(od(j).LFP,2);
	Shat = filter(ones(300,1),300,S);

	[stim_peaks, stim_peak_loc] = findpeaks(S,'MinPeakProminence',.25,'MinPeakDistance',100);
	min_stim_before_whiff = NaN*stim_peak_loc;
	min_stim_before_whiff_locs = NaN*stim_peak_loc;
	peak_firing_rate = NaN*stim_peaks;
	peak_LFP = NaN*stim_peaks;
	peak_firing_locs = NaN*stim_peaks;
	peak_LFP_locs = NaN*stim_peaks;

	for i = 2:(length(stim_peak_loc)-1)
		a = max([stim_peak_loc(i) - 90 stim_peak_loc(i-1)]);
		[min_stim_before_whiff(i), loc] = min(S(a:stim_peak_loc(i)));
		min_stim_before_whiff_locs(i) = loc + a;
	end
	for i = 2:(length(stim_peak_loc)-1)
		z = min([stim_peak_loc(i)+200 min_stim_before_whiff_locs(i+1)]);
		[peak_firing_rate(i),loc] = max(R(stim_peak_loc(i):z));
		peak_firing_locs(i) = loc + stim_peak_loc(i);

		z = min([stim_peak_loc(i)+200 min_stim_before_whiff_locs(i+1)]);
		[peak_LFP(i),loc] = min(X(stim_peak_loc(i):z));
		peak_LFP_locs(i) = loc + stim_peak_loc(i);
	end

	firing_rate_before_whiff = [NaN; R(nonnans(min_stim_before_whiff_locs)); NaN];
	LFP_before_whiff = [NaN; X(nonnans(min_stim_before_whiff_locs)); NaN];
	stim_300ms_before_whiff = [NaN; Shat(nonnans(min_stim_before_whiff_locs)); NaN];

	orn = j*ones(length(LFP_before_whiff),1);

	% assemble
	f = fieldnames(whiffs);
	for i = 1:length(f)
		eval(['whiffs.' f{i} '= [whiffs.' f{i} '(:);  '  f{i} '(:)];']);
	end

end

% check
figure('outerposition',[0 0 1300 902],'PaperUnits','points','PaperSize',[1300 902]); hold on
subplot(3,4,1); hold on
c = lines(5);
for i = 1:5
	plot(whiffs.stim_peaks(whiffs.orn==i),whiffs.peak_LFP(whiffs.orn==i),'+','Color',c(i,:))
end
xlabel('Peak stimulus (V)')
ylabel('Peak LFP response (mV)')

% same plot, but colour by LFP before whiff 
subplot(3,4,2); hold on
cinfo = whiffs.LFP_before_whiff;
cinfo = cinfo - min(cinfo);
cinfo = cinfo/max(cinfo);
cinfo = floor(1+cinfo*99);
ignore_these = isnan(cinfo);
cinfo(ignore_these) = 1;
c = parula(100);
cc = c(cinfo,:);
scatter(whiffs.stim_peaks,whiffs.peak_LFP,24,cc,'filled')
ch = colorbar('east');
ch.Position = [0.4700    0.79    0.0154    0.13];
caxis([min(whiffs.LFP_before_whiff) max(whiffs.LFP_before_whiff)])
xlabel('Peak stimulus (V)')
ylabel('Peak LFP response (mV)')
title('Color is LFP before whiff (mV)')

subplot(3,4,3); hold on
c = lines(5);
for i = 1:5
	plot(whiffs.stim_peaks(whiffs.orn==i),whiffs.peak_firing_rate(whiffs.orn==i),'+','Color',c(i,:))
end
xlabel('Peak stimulus (V)')
ylabel('Peak firing rate (Hz)')

% same plot, but colour by LFP before whiff 
subplot(3,4,4); hold on
cinfo = whiffs.LFP_before_whiff;
cinfo = cinfo - min(cinfo);
cinfo = cinfo/max(cinfo);
cinfo = floor(1+cinfo*99);
ignore_these = isnan(cinfo);
cinfo(ignore_these) = 1;
c = parula(100);
cc = c(cinfo,:);
scatter(whiffs.stim_peaks,whiffs.peak_firing_rate,24,cc,'filled')
ch = colorbar('east');
ch.Position = [0.8772    0.7288    0.0140    0.0738];
caxis([min(whiffs.LFP_before_whiff) max(whiffs.LFP_before_whiff)])
xlabel('Peak stimulus (V)')
ylabel('Peak firing rate (Hz)')
title('Color is LFP before whiff (mV)')


subplot(3,4,5); hold on
c = lines(5);
for i = 1:5
	plot(whiffs.LFP_before_whiff(whiffs.orn==i),whiffs.peak_firing_rate(whiffs.orn==i),'+','Color',c(i,:))
end
xlabel('LFP before whiff (mV)')
ylabel('Peak firing rate (Hz)')

subplot(3,4,6); hold on
c = lines(5);
for i = 1:5
	x = whiffs.LFP_before_whiff(whiffs.orn==i);
	y = whiffs.peak_firing_rate(whiffs.orn==i) -  whiffs.firing_rate_before_whiff(whiffs.orn==i);
	plot(x,y,'+','Color',c(i,:))
end
xlabel('LFP before whiff (mV)')
ylabel('\Delta firing rate (Hz)')

subplot(3,4,7); hold on
c = lines(5);
for i = 1:5
	x = whiffs.LFP_before_whiff(whiffs.orn==i);
	y = whiffs.peak_firing_rate(whiffs.orn==i) -  whiffs.firing_rate_before_whiff(whiffs.orn==i);
	y = y./(whiffs.stim_peaks(whiffs.orn==i) - whiffs.min_stim_before_whiff(whiffs.orn==i));
	plot(x,y,'+','Color',c(i,:))
end
xlabel('LFP before whiff (mV)')
ylabel('\DeltaR/\DeltaS (Hz/V)')

subplot(3,4,8); hold on
c = lines(5);
for i = 1:5
	x = whiffs.LFP_before_whiff(whiffs.orn==i);
	y = whiffs.peak_firing_rate(whiffs.orn==i);
	y = y./whiffs.stim_peaks(whiffs.orn==i);
	plot(x,y,'+','Color',c(i,:))
end
xlabel('LFP before whiff (mV)')
ylabel('R/S (Hz/V)')

subplot(3,4,9); hold on
c = lines(5);
for i = 1:5
	x = whiffs.min_stim_before_whiff(whiffs.orn==i);
	y = whiffs.peak_firing_rate(whiffs.orn==i) -  whiffs.firing_rate_before_whiff(whiffs.orn==i);
	y = y./(whiffs.stim_peaks(whiffs.orn==i) - whiffs.min_stim_before_whiff(whiffs.orn==i));
	plot(x,y,'+','Color',c(i,:))
end
xlabel('min S before whiff (V)')
ylabel('\DeltaR/\DeltaS (Hz/V)')
set(gca,'XScale','log','YScale','log','XTick',[1e-3 1e-2 1e-1 1 1e1],'YTick',[1e-1 1 10 100 1e3],'XLim',[1e-3 10],'YLim',[1e-1 1e3])

subplot(3,4,10); hold on
c = lines(5);
for i = 1:5
	x = whiffs.stim_300ms_before_whiff(whiffs.orn==i);
	y = whiffs.peak_firing_rate(whiffs.orn==i) -  whiffs.firing_rate_before_whiff(whiffs.orn==i);
	y = y./(whiffs.stim_peaks(whiffs.orn==i) - whiffs.min_stim_before_whiff(whiffs.orn==i));
	plot(x,y,'+','Color',c(i,:))
end
xlabel('<S> 300 ms before whiff (V)')
ylabel('\DeltaR/\DeltaS (Hz/V)')
set(gca,'XScale','log','YScale','log','XTick',[1e-3 1e-2 1e-1 1 1e1],'YTick',[1e-1 1 10 100 1e3],'XLim',[1e-3 10],'YLim',[1e-1 1e3])

subplot(3,4,11); hold on
c = lines(5);
for i = 1:5
	x = whiffs.stim_300ms_before_whiff(whiffs.orn==i);
	y = whiffs.peak_LFP(whiffs.orn==i) -  whiffs.LFP_before_whiff(whiffs.orn==i);
	y = y./(whiffs.stim_peaks(whiffs.orn==i) - whiffs.min_stim_before_whiff(whiffs.orn==i));
	y = abs(y);
	plot(x,y,'+','Color',c(i,:))
end
xlabel('<S> 300 ms before whiff (V)')
ylabel('abs(\DeltaLFP)/\DeltaS (Hz/V)')
set(gca,'XScale','log','YScale','log','XTick',[1e-3 1e-2 1e-1 1 1e1],'YTick',[1e-2 1e-1 1 10 100 1e3],'XLim',[1e-3 10],'YLim',[1e-2 1e2])

subplot(3,4,12); hold on
c = lines(5);
for i = 1:5
	x = whiffs.min_stim_before_whiff(whiffs.orn==i);
	y = whiffs.peak_firing_rate(whiffs.orn==i);
	y = y./whiffs.stim_peaks(whiffs.orn==i);
	plot(x,y,'+','Color',c(i,:))
end
xlabel('min S before whiff (V)')
ylabel('R/S (Hz/V)')
set(gca,'XScale','log','YScale','log','XTick',[1e-3 1e-2 1e-1 1 1e1],'YTick',[1e-1 1 10 100 1e3],'XLim',[1e-3 10],'YLim',[1e-1 1e3])

prettyFig();

if being_published
	snapnow
	delete(gcf)
end




