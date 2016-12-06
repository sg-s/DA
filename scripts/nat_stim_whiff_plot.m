


pHeader;





%% Detailed look at naturalistic stimulus whiffs
% In this document, I look at the ORN responses to individual whiffs in the naturalistic stimulus response, to try to understand what's going on in the response without any filter analysis. 


% get the naturalistic stimuli 
clear ab3 ab2
load(getPath(dataManager,'5c7dacc5b42ff0eebb980d80fec120c3'),'data','spikes')
PID = data(2).PID;
time = 1e-4*(1:length(PID));
all_spikes = spikes(2).A;

% A spikes --> firing rate
fA = spiketimes2f(all_spikes,time);

tA = 1e-3*(1:length(fA));
PID2 = fA;
for i = 1:width(PID2)
	PID2(:,i) = interp1(time,PID(i,:),tA);
end
PID = PID2; clear PID2
% some minor cleaning up
PID(end,:) = PID(end-1,:); 

% remove the baseline from the PID, and remember the error
PID_baseline = mean(mean(PID(1:5e3,:)));
PID = PID - PID_baseline;

NSdata.fA = mean(fA,2);
NSdata.PID = mean(PID,2);
NSdata.time = 1e-3*(1:length(NSdata.fA));

S = NSdata.PID;
R = NSdata.fA;

figure('outerposition',[0 0 1541 902],'PaperUnits','points','PaperSize',[1541 902]); hold on

a = [8.901 6.124; 26.09 58.04];
before = .5;
after = .5;

c = lines(2);

% show the stimulus
subplot(3,4,1:4); hold on
plot(tA,S,'k')
ylabel('Stimulus (V)')
set(gca,'XLim',[0 70],'YLim',[-.2 9])

% indicate the segments we zoom into we use
plot([a(1,1)-before a(1,1)+after],[4 4],'Color',c(1,:),'LineWidth',4)
plot([a(1,2)-before a(1,2)+after],[4 4],'Color',c(2,:),'LineWidth',4)

plot([a(2,1)-before a(2,1)+after],[8 8],'Color',c(1,:),'LineWidth',4)
plot([a(2,2)-before a(2,2)+after],[8 8],'Color',c(2,:),'LineWidth',4)

% show response
subplot(3,4,5:8); hold on
plot(tA,R,'k')
ylabel('ab3A firing rate (Hz)')
set(gca,'XLim',[0 70],'YLim',[-5 130])

% indicate the segments we zoom into we use
plot([a(1,1)-before a(1,1)+after],[125 125],'Color',c(1,:),'LineWidth',4)
plot([a(1,2)-before a(1,2)+after],[125 125],'Color',c(2,:),'LineWidth',4)

plot([a(2,1)-before a(2,1)+after],[125 125],'Color',c(1,:),'LineWidth',4)
plot([a(2,2)-before a(2,2)+after],[125 125],'Color',c(2,:),'LineWidth',4)


for i = 1:2

	for j = 1:2
		subplot(3,4,8+(i-1)*2+1); hold on
		s = 1e3*(a(i,j) -  before);
		e = 1e3*(a(i,j) + after);
		x = tA(s:e); x = x - a(i,j);
		y = S(s:e);

		plot(x,y,'Color',c(j,:))

		subplot(3,4,8+(i-1)*2+2); hold on
		s = 1e3*(a(i,j) -  before);
		e = 1e3*(a(i,j) + after);
		x = tA(s:e); x = x - a(i,j);
		y = R(s:e);

		plot(x,y,'Color',c(j,:))

	end
end

subplot(3,4,9); hold on
xlabel('Time (s)')
ylabel('Stimulus (V)')

subplot(3,4,11); hold on
ylabel('Stimulus (V)')

subplot(3,4,10); hold on
ylabel('Firing rate (Hz) (V)')

subplot(3,4,12); hold on
ylabel('Firing rate (Hz) (V)')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%% Version Info
%
pFooter;

