function [] = plotScaledNatStimData(data,being_published,t)
time = (1:length(data(1).S))*1e-3;
for i = 1:length(data)
	figure('outerposition',[0 0 1000 901],'PaperUnits','points','PaperSize',[1000 901]); hold on
	ax(1) = subplot(3,1,1); hold on
	ax(2) = subplot(3,1,2); hold on
	ax(3) = subplot(3,1,3); hold on
	c = lines(size(data(i).S,2));
	for j = 1:size(data(i).S,2)
		plot(ax(1),time,data(i).S(:,j),'Color',c(j,:))
		plot(ax(2),time,data(i).X(:,j),'Color',c(j,:))
		plot(ax(3),time,data(i).R(:,j),'Color',c(j,:))
	end

	xlabel(ax(3),'Time (s)')
	ylabel(ax(1),'Stimulus (V)')
	ylabel(ax(2),'\DeltaLFP (mV)')
	ylabel(ax(3),'Firing rate (Hz)')

	suptitle(t)

	prettyFig();
	
	if being_published
		snapnow
		delete(gcf)
	end
end