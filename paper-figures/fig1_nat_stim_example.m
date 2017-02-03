% fig1_nat_stim_example

%%
% This script makes fig 1 for the eLife resubmission. The point of this figure is to show our naturalistic stimulus data, and the features of the ORN response we see in it. 

pHeader;


% make the figure and all subplots 
figure('outerposition',[0 0 1302 1001],'PaperUnits','points','PaperSize',[1000 1001]); hold on
clear ax
ax(1) = subplot(4,4,1:3); hold on
ax(2) = subplot(4,4,5:7); hold on
ax(3) = subplot(4,4,9:11); hold on

ax(4) = subplot(4,4,13); hold on
ax(5) = subplot(4,4,14); hold on

ax(6) = subplot(4,4,4); hold on
ax(7) = subplot(4,4,8); hold on
ax(8) = subplot(4,4,12); hold on


ax(9) = subplot(4,4,15); hold on
ax(10) = subplot(4,4,16); hold on


% get the data

cdata = consolidateData2(getPath(dataManager,'4608c42b12191b383c84fea52392ea97'));
[cdata, data] =  assembleScaledNatStim(cdata);
time = 1e-3*(1:length(data(1).S));

% plot the data
c = parula(4); 
for i = 1:3
	plot(ax(1),time,data(2).S(:,i),'Color',c(i,:))
	plot(ax(2),time,data(2).X(:,i),'Color',c(i,:))
	plot(ax(3),time,data(2).R(:,i),'Color',c(i,:))
end

% show whiff statistics 
i = 2;
for j = 1:size(data(i).S,2)
	S = data(i).S(:,j);
	X = data(i).X(:,j);
	R = data(i).R(:,j);
	ws = whiffStatistics(S,X,R,300,'MinPeakProminence',max(S/1e2),'debug',false);
	plot(ax(4),ws.stim_peaks,ws.peak_LFP,'.','MarkerSize',20,'Color',c(j,:))
	plot(ax(5),ws.stim_peaks,ws.peak_firing_rate,'.','MarkerSize',20,'Color',c(j,:))
end

% show context-dependent variation in whiff response
s_range = [.8 1.2];

show_these = [2       27764
           2       28790
           2       59776
           3       58049];

for i = 2
	for j = 1:length(show_these)
		this_stim = show_these(j,1);
		this_loc = show_these(j,2);

		S = data(i).S(:,this_stim);
		X = data(i).X(:,this_stim);
		R = data(i).R(:,this_stim);

		a = this_loc - 300;
		z = this_loc+300;

		t = (1:length(S(a:z))) - 300;

		plot(ax(6),t,S(a:z))
		plot(ax(7),t,X(a:z))
		plot(ax(8),t,R(a:z))

	end
end

% show linear analysis output plot


ss = 50;

i = 2;
for j = 1:size(data(i).S,2)
	S = data(i).S(:,j);
	X = data(i).X(:,j);
	R = data(i).R(:,j);

	K1 = fitFilter2Data(S,X,'offset',200);
	K2 = fitFilter2Data(S,R,'offset',200);

	K1 = K1(100:end-100);
	K2 = K2(100:end-100);
	filtertime = 1e-3*(1:length(K1)) - .1;
	Xp = convolve(time,S,K1,filtertime);
	Rp = convolve(time,S,K2,filtertime);
	K1 = K1/(nanstd(Xp)/nanstd(S)); % normalise correctly 
	Xp = convolve(time,S,K1,filtertime);
	K2 = K2/(nanstd(Rp)/nanstd(S)); % normalise correctly 
	Rp = convolve(time,S,K2,filtertime);

 
	plot(ax(9),Xp(1:ss:end),X(1:ss:end),'.-','Color',c(:,j))
	plot(ax(10),Rp(1:ss:end),R(1:ss:end),'.-','Color',c(:,j))


end


% axes modifications, labels, etc.
ylabel(ax(1),'Stimulus (V)')
set(ax(1),'XLim',[0 70])
ylabel(ax(2),'\Delta LFP (mV)')
set(ax(2),'XLim',[0 70])
ylabel(ax(3),['ab2A firing' char(10) 'rate (Hz)'])
xlabel(ax(3),'Time (s)')
set(ax(3),'XLim',[0 70])

xlabel(ax(4),'Whiff intensity (V)')
ylabel(ax(4),'Peak LFP (mV)')
set(ax(4),'XScale','log','XLim',[1e-2 1e1],'XTick',[1e-2 1e-1 1e0 1e1])
axes(ax(4))
axis square

xlabel(ax(5),'Whiff intensity (V)')
ylabel(ax(5),'Peak firing rate (Hz)')
set(ax(5),'XScale','log','YLim',[0 300],'XLim',[1e-2 1e1],'XTick',[1e-2 1e-1 1e0 1e1])
axes(ax(5))
axis square

xlabel(ax(8),'Time since whiff (ms)')

xlabel(ax(9),'Projected stimulus (V)')
ylabel(ax(9),'LFP (mV)')
set(ax(9),'YDir','reverse','XDir','reverse')
axes(ax(9))
axis square

xlabel(ax(10),'Projected stimulus (V)')
ylabel(ax(10),'Firing rate (Hz)')
axes(ax(10))
axis square


% move some things around
ax(1).Position = [0.1300    0.8    0.5689    0.11];
ax(2).Position = [ 0.1300    0.625    0.5689    0.11];
ax(3).Position = [ 0.1300    0.475    0.5689    0.11];
for i = 6:8
	ax(i).Position(2) = ax(i-5).Position(2);
	ax(i).Position(4) = ax(i-5).Position(4);
end

for i = [4 5 9 10]
	movePlot(ax(i),'left',.04)
	ax(i).Position(2) = .15;
	ax(i).Position(3:4) = .22;
end



prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end


%% Version Info
%
pFooter;


