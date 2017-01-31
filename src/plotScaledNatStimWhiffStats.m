function  [] = plotScaledNatStimWhiffStats(data)

w = min([1300 500*length(data)]);

figure('outerposition',[0 0 w 901],'PaperUnits','points','PaperSize',[1300 901]); hold on

for i = 1:length(data)
	c = lines(size(data(i).S,2));
	for j = 1:size(data(i).S,2)
		S = data(i).S(:,j);
		X = data(i).X(:,j);
		R = data(i).R(:,j);
		ws = whiffStatistics(S,X,R,300,'MinPeakProminence',max(S/1e2),'debug',false);

		subplot(2,length(data),i); hold on
		plot(ws.stim_peaks,ws.peak_LFP,'.','MarkerSize',20,'Color',c(j,:))
		set(gca,'XScale','log')
		ylabel('LFP_{peak} (mV)')
		axis square

		subplot(2,length(data),length(data)+i); hold on
		plot(ws.stim_peaks,ws.peak_firing_rate,'.','MarkerSize',20,'Color',c(j,:))
		set(gca,'XScale','log','YLim',[0 300])
		xlabel('Stim_{peak} (V)')
		ylabel('Firing rate_{peak} (mV)')
		axis square
	end
end