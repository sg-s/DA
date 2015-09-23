function [] = PlotFilterDiagnostics2(diagnostics,marker_size,marker_size2,font_size,side_string)
bf = diagnostics.bestfilter;
figure('outerposition',[0 0 1200 400],'PaperUnits','points','PaperSize',[1200 400]); hold on
subplot(1,5,1), hold on
plot(diagnostics.reg,diagnostics.err,'.','MarkerSize',marker_size2), hold on
plot(diagnostics.reg(bf),diagnostics.err(bf),'r.','MarkerSize',marker_size2)
set(gca,'LineWidth',2,'FontSize',20,'box','on','XScale','log')
xlabel('r','FontSize',font_size)
ylabel(side_string,'FontSize',font_size)
title('r-square','FontSize',font_size)
axis square

subplot(1,5,2), hold on
plot(diagnostics.reg,diagnostics.filter_sum,'.','MarkerSize',marker_size2)
plot(diagnostics.reg(bf),diagnostics.filter_sum(bf),'r.','MarkerSize',marker_size2)
set(gca,'LineWidth',2,'FontSize',20,'box','on','XScale','log')
xlabel('r','FontSize',font_size)
title('\Sigma|diff(K)|','FontSize',font_size)
if max(diagnostics.filter_sum)/min(diagnostics.filter_sum) > 100
	% use a log axis
	set(gca,'YScale','log')
end
axis square

subplot(1,5,3), hold on
plot(diagnostics.reg,diagnostics.filter_height,'.','MarkerSize',marker_size2)
plot(diagnostics.reg(bf),diagnostics.filter_height(bf),'r.','MarkerSize',marker_size2)
set(gca,'LineWidth',2,'FontSize',20,'box','on','XScale','log')
xlabel('r','FontSize',font_size)
title('Filter Height','FontSize',font_size)
if max(diagnostics.filter_height)/min(diagnostics.filter_height) > 100
	% use a log axis
	set(gca,'YScale','log')
end
axis square

subplot(1,5,4), hold on
plot(diagnostics.reg,diagnostics.slope,'.','MarkerSize',marker_size2)
plot(diagnostics.reg(bf),diagnostics.slope(bf),'r.','MarkerSize',marker_size2)
set(gca,'LineWidth',2,'FontSize',20,'box','on','XScale','log')
xlabel('r','FontSize',font_size)
title('Gain','FontSize',font_size)
if max(diagnostics.slope)/min(diagnostics.slope) > 100
	% use a log axis
	set(gca,'YScale','log')
end
axis square

subplot(1,5,5), hold on
plot(diagnostics.reg,diagnostics.mcond,'.','MarkerSize',marker_size2)
plot(diagnostics.reg(bf),diagnostics.mcond(bf),'r.','MarkerSize',marker_size2)
set(gca,'LineWidth',2,'FontSize',20,'box','on','XScale','log')
xlabel('r','FontSize',font_size)
title('C#','FontSize',font_size)
if max(diagnostics.mcond)/min(diagnostics.mcond) > 100
	% use a log axis
	set(gca,'YScale','log')
end
axis square
