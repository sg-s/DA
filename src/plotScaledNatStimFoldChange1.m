% helper function for nat_stim_analysis

function [fold_changes] = plotScaledNatStimFoldChange1(data,being_published)

clear fold_changes

figure('outerposition',[0 0 1000 901],'PaperUnits','points','PaperSize',[1000 901]); hold on
for i = 1:length(data)
	S0 = data(i).S(:,end);
	X0 = data(i).X(:,end);
	R0 = data(i).R(:,end);
	c = lines(size(data(i).S,2));
	for j = 1:size(data(i).S,2)-1
		subplot(3,length(data),i); hold on
		xlabel('Stimulus')
		ylabel('Stimulus')
		title(['ORN #' oval(i)])
		plot(S0,data(i).S(:,j),'Color',c(j,:))
		% compute stimulus fold change
		ff = fit(S0,data(i).S(:,j),'poly1');
		fold_changes(i).S(j) = ff.p1;
		cf = confint(ff);
		fold_changes(i).S_l(j) = cf(1,1);
		fold_changes(i).S_h(j) = cf(2,1);
	end
	set(gca,'XLim',[0 1.5*max(max([data.S]))],'YLim',[0 1.5*max(max([data.S]))])
	plot([0 1.5*max(max([data.S]))], [0  1.5*max(max([data.S]))],'k--')
	axis square

	for j = 1:size(data(i).S,2)-1
		subplot(3,length(data),length(data)+i); hold on
		plot(X0,data(i).X(:,j),'Color',c(j,:))
		xlabel('LFP')
		ylabel('LFP')
		% compute LFP fold change
		ff = fit(X0,data(i).X(:,j),'poly1');
		fold_changes(i).X(j) = ff.p1;
		cf = confint(ff);
		fold_changes(i).X_l(j) = cf(1,1) - ff.p1;
		fold_changes(i).X_h(j) = ff.p1 - cf(2,1);
	end
	%set(gca,'XLim',[0 1.5*max(max([data.S]))],'YLim',[0 1.5*max(max([data.S]))])
	axis square

	for j = 1:size(data(i).S,2)-1
		subplot(3,length(data),2*length(data)+i); hold on
		plot(R0,data(i).R(:,j),'Color',c(j,:))
		xlabel('Firing rate')
		ylabel('Firing rate')
		fold_changes(i).R = 0*fold_changes(i).X;
		fold_changes(i).R_l = 0*fold_changes(i).X;
		fold_changes(i).R_h = 0*fold_changes(i).X;
		% compute firing fold change
		try
			ff = fit(R0,data(i).R(:,j),'poly1');
			fold_changes(i).R(j) = ff.p1;
			cf = confint(ff);
			fold_changes(i).R_l(j) = cf(1,1) - ff.p1;
			fold_changes(i).R_h(j) = ff.p1 - cf(2,1);
		catch
		end
	end
	%set(gca,'XLim',[0 1.5*max(max([data.S]))],'YLim',[0 1.5*max(max([data.S]))])
	axis square

end
prettyFig();

if being_published
	snapnow
	delete(gcf)
end