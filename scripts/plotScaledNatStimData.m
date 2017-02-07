function [] = plotScaledNatStimData(data,being_published,t,make_plot)

ss = 100;

if nargin < 4
	make_plot = true;
end

time = (1:length(data(1).S))*1e-3;
for i = 1:length(data)
	if make_plot
		figure('outerposition',[0 0 1000 901],'PaperUnits','points','PaperSize',[1000 901]); hold on
		ax(1) = subplot(3,5,1:5); hold on
		ax(2) = subplot(3,5,6:10); hold on
		ax(3) = subplot(3,5,11:15); hold on
	end

	% ax_dist(1) = subplot(3,5,5); hold on; 
	% ax_dist(2) = subplot(3,5,10); hold on
	% ax_dist(3) = subplot(3,5,15); hold on
	% for j = 1:3
	% 	set(ax_dist(j),'YScale','log','XScale','log')
	% end
	c = lines(size(data(i).S,2));
	for j = 1:size(data(i).S,2)
		if make_plot
			plot(ax(1),time(1:ss:end),data(i).S(1:ss:end,j),'Color',c(j,:))
			plot(ax(2),time(1:ss:end),data(i).X(1:ss:end,j),'Color',c(j,:))
			plot(ax(3),time(1:ss:end),data(i).R(1:ss:end,j),'Color',c(j,:))
		end
	end

	if make_plot
		xlabel(ax(3),'Time (s)')
		ylabel(ax(1),'Stimulus (V)')
		ylabel(ax(2),'\DeltaLFP (mV)')
		ylabel(ax(3),'Firing rate (Hz)')
		set(ax(1),'XLim',[0 70])
		set(ax(2),'XLim',[0 70])
		set(ax(3),'XLim',[0 70])

		suptitle(t)

		prettyFig();
	
		if being_published
			snapnow
			delete(gcf)
		end
	end
end