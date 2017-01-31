function [] = plotScaledNatStimFoldChange2(fold_changes,being_published)


figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
c = lines(length(fold_changes));
clear l
for i = 1:length(fold_changes)
	plot(fold_changes(i).S,fold_changes(i).X,'LineStyle','none','Color',c(i,:),'Marker','.','MarkerSize',20)
	plot(fold_changes(i).S,fold_changes(i).R,'LineStyle','none','Color',c(i,:),'Marker','d','MarkerSize',10)
end
plot([0 10],[0 10],'k--')
l(1) = plot(NaN,NaN,'.','MarkerSize',20,'Color','k');
l(2) = plot(NaN,NaN,'d','MarkerSize',10,'Color','k');
legend(l,{'LFP','Firing rate'})
set(gca,'YLim',[.5 3])
xlabel('Stimulus scale')
ylabel('Response scale')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end