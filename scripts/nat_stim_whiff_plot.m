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

figure('outerposition',[0 0 1541 902],'PaperUnits','points','PaperSize',[1541 902]); hold on

a = [36.36 47.9; 47.9 59.46];
before = .5;
after = .5;

c = lines(2);

% show the stimulus
subplot(3,5,1:3); hold on
plot(tA,S,'k')
ylabel('Stimulus (V)')
set(gca,'XLim',[0 70],'YLim',[-.2 9])

% indicate the segments we zoom into we use
plot([a(1,1)-before a(1,1)+after],[2 2],'Color',c(1,:),'LineWidth',4)
plot([a(1,2)-before a(1,2)+after],[2 2],'Color',c(2,:),'LineWidth',4)

plot([a(2,1)-before a(2,1)+after],[1.5 1.5],'Color',c(1,:),'LineWidth',4)
plot([a(2,2)-before a(2,2)+after],[1.5 1.5],'Color',c(2,:),'LineWidth',4)

% show LFP response
subplot(3,5,6:8); hold on
plot(tA,X,'k')
ylabel('ab3 LFP (\DeltamV)')
set(gca,'XLim',[0 70],'YLim',[-20 5])

% % indicate the segments we zoom into we use
% plot([a(1,1)-before a(1,1)+after],[125 125],'Color',c(1,:),'LineWidth',4)
% plot([a(1,2)-before a(1,2)+after],[125 125],'Color',c(2,:),'LineWidth',4)

% plot([a(2,1)-before a(2,1)+after],[125 125],'Color',c(1,:),'LineWidth',4)
% plot([a(2,2)-before a(2,2)+after],[125 125],'Color',c(2,:),'LineWidth',4)

% show the firing rate
subplot(3,5,11:13); hold on
plot(tA,R,'k')
ylabel('ab3A firing rate (Hz)')
set(gca,'XLim',[0 70],'YLim',[-5 170])

% % indicate the segments we zoom into we use
% plot([a(1,1)-before a(1,1)+after],[125 125],'Color',c(1,:),'LineWidth',4)
% plot([a(1,2)-before a(1,2)+after],[125 125],'Color',c(2,:),'LineWidth',4)

% plot([a(2,1)-before a(2,1)+after],[125 125],'Color',c(1,:),'LineWidth',4)
% plot([a(2,2)-before a(2,2)+after],[125 125],'Color',c(2,:),'LineWidth',4)


stim_plots = [4 5];
lfp_plots = [9 10];
firing_plots = [14 15];
for i = 1:2
	for j = 1:2
		% plot the stimulus
		subplot(3,5,stim_plots(i)); hold on
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
		subplot(3,5,lfp_plots(i)); hold on
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
		subplot(3,5,firing_plots(i)); hold on
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

subplot(3,5,4); hold on
set(gca,'YLim',[-.1 .7])

subplot(3,5,5); hold on
set(gca,'YLim',[-.1 .7])

subplot(3,5,11:13); hold on
xlabel('Time (s)')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%% 
% can we also show the slowdown in the LFP and the constancy of the firing rate response? 


clear od
load('/local-data/DA-paper/data-for-paper/fig7/nat-stim-ab3/combined_data.ORNData','-mat')
DNSdata.LFP = [];
DNSdata.fA = [];
DNSdata.PID = [];
DNSdata.orn = [];
for i = length(od):-1:1
	DNSdata.LFP = [DNSdata.LFP  od(i).LFP];
	DNSdata.orn = [DNSdata.orn; vectorise(i*ones(od(i).n_trials,1))];
	DNSdata.fA = [DNSdata.fA  od(i).firing_rate];
	DNSdata.PID = [DNSdata.PID  od(i).stimulus];
end
% remove baseline from stimulus
for i = 1:size(DNSdata.PID,2)
	DNSdata.PID(:,i) = DNSdata.PID(:,i) - min(DNSdata.PID(1:5e3,i));
end


S = nanmean(DNSdata.PID(:,DNSdata.orn == 6),2);
R = nanmean(DNSdata.fA(:,DNSdata.orn == 6),2);
X = nanmean(DNSdata.LFP(:,DNSdata.orn == 6),2);

% find all peaks in the stimulus
[stim_peaks,loc] = findpeaks(S,'MinPeakProminence',.25,'MinPeakDistance',200);

% for each, find a peak in the LFP corresponding to these peaks
xloc = NaN*loc;
for i = 1:length(loc)
	[~,xloc(i)]= max(-X(loc(i):loc(i)+500));
	xloc(i) = xloc(i) + loc(i);
end

% also find a peak in the firing rate
rloc = NaN*loc;
for i = 1:length(loc)
	[~,rloc(i)]= max(R(loc(i):loc(i)+500));
	rloc(i) = rloc(i) + loc(i);
end

figure('outerposition',[0 0 1501 500],'PaperUnits','points','PaperSize',[1501 500]); hold on
subplot(1,3,1); hold on
plot_this = [ 12     9    46    36    27];
set(gca,'ColorOrder',parula(length(plot_this)));
peak_stim = [];
for i = plot_this
	peak_stim = [peak_stim max(S(loc(i)-300:loc(i)+300))];
	plot(S(loc(i)-300:loc(i)+300))
end


subplot(1,3,3); hold on
set(gca,'ColorOrder',parula(length(plot_this)));
peak_response = [];
for i = plot_this
	temp = R(loc(i)-300:loc(i)+300);
	peak_response = [peak_response max(temp)];
	temp = temp - mean(temp(1:200));
	plot(temp/max(temp))

end
set(gca,'YLim',[0 1])

subplot(1,3,2); hold on
set(gca,'ColorOrder',parula(length(plot_this)));
for i = plot_this
	temp = X(loc(i)-300:loc(i)+400);
	temp = temp - mean(temp(1:200));
	temp = -temp;
	plot(temp/max(temp))
end
set(gca,'YLim',[0 1])

prettyFig();

if being_published
	snapnow
	delete(gcf)
end


%% Version Info
%
pFooter;

