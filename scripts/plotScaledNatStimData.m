function [] = plotScaledNatStimData(data,being_published,t)

ss = 100;

time = (1:length(data(1).S))*1e-3;
for i = 1:length(data)
	figure('outerposition',[0 0 1000 901],'PaperUnits','points','PaperSize',[1000 901]); hold on
	ax(1) = subplot(3,5,1:5); hold on
	ax(2) = subplot(3,5,6:10); hold on
	ax(3) = subplot(3,5,11:15); hold on

	% ax_dist(1) = subplot(3,5,5); hold on; 
	% ax_dist(2) = subplot(3,5,10); hold on
	% ax_dist(3) = subplot(3,5,15); hold on
	% for j = 1:3
	% 	set(ax_dist(j),'YScale','log','XScale','log')
	% end
	c = lines(size(data(i).S,2));
	for j = 1:size(data(i).S,2)
		plot(ax(1),time(1:ss:end),data(i).S(1:ss:end,j),'Color',c(j,:))

		% also plot the distribution 
		% temp = data(i).S(:,j);
		% [hy,hx] = histcounts(temp,100); hy = hy/sum(hy);
		% hx = hx(1:end-1) + mean(diff(hx));
		% plot(ax_dist(1),hx,hy,'+','Color',c(j,:))

		plot(ax(2),time(1:ss:end),data(i).X(1:ss:end,j),'Color',c(j,:))

		% temp = -data(i).X(:,j);
		% [hy,hx] = histcounts(temp,100); hy = hy/sum(hy);
		% hx = hx(1:end-1) + mean(diff(hx));
		% plot(ax_dist(2),hx,hy,'+','Color',c(j,:))

		plot(ax(3),time(1:ss:end),data(i).R(1:ss:end,j),'Color',c(j,:))

		% temp = data(i).R(:,j);
		% [hy,hx] = histcounts(temp,100); hy = hy/sum(hy);
		% hx = hx(1:end-1) + mean(diff(hx));
		% plot(ax_dist(3),hx,hy,'+','Color',c(j,:))
	end

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