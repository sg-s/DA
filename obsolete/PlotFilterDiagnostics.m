% PlotFilterDiagnostics
% plots data from the diagnostics structure created by FitFitler2Data
% 
% created by Srinivas Gorur-Shandilya at 2:50 , 10 July 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [] = PlotFilterDiagnostics(diagnostics,marker_size,marker_size2,font_size)
bf = diagnostics.C.bestfilter;
figure('outerposition',[0 0 1200 800],'PaperUnits','points','PaperSize',[1200 800]); hold on
subplot(2,4,1), hold on
plot(diagnostics.C.reg,diagnostics.C.err,'.','MarkerSize',marker_size2), hold on
plot(diagnostics.C.reg(bf),diagnostics.C.err(bf),'r.','MarkerSize',marker_size2)
set(gca,'LineWidth',2,'FontSize',20,'box','on','XScale','log')
xlabel('r','FontSize',font_size)
ylabel('Chichilnisky','FontSize',font_size)
title('Error','FontSize',font_size)
axis square

subplot(2,4,2), hold on
plot(diagnostics.C.reg,diagnostics.C.filter_sum,'.','MarkerSize',marker_size2)
plot(diagnostics.C.reg(bf),diagnostics.C.filter_sum(bf),'r.','MarkerSize',marker_size2)
set(gca,'LineWidth',2,'FontSize',20,'box','on','XScale','log')
xlabel('r','FontSize',font_size)
title('\Sigma|K|','FontSize',font_size)
axis square

subplot(2,4,3), hold on
plot(diagnostics.C.reg,diagnostics.C.filter_height,'.','MarkerSize',marker_size2)
plot(diagnostics.C.reg(bf),diagnostics.C.filter_height(bf),'r.','MarkerSize',marker_size2)
set(gca,'LineWidth',2,'FontSize',20,'box','on','XScale','log')
xlabel('r','FontSize',font_size)
title('Filter Height','FontSize',font_size)
axis square

subplot(2,4,4), hold on
plot(diagnostics.C.reg,diagnostics.C.slope,'.','MarkerSize',marker_size2)
plot(diagnostics.C.reg(bf),diagnostics.C.slope(bf),'r.','MarkerSize',marker_size2)
set(gca,'LineWidth',2,'FontSize',20,'box','on','XScale','log')
xlabel('r','FontSize',font_size)
title('Gain','FontSize',font_size)
axis square

bf = diagnostics.D.bestfilter;
subplot(2,4,5), hold on
plot(diagnostics.D.reg,diagnostics.D.err,'.','MarkerSize',marker_size2)
plot(diagnostics.D.reg(bf),diagnostics.D.err(bf),'r.','MarkerSize',marker_size2)
set(gca,'LineWidth',2,'FontSize',20,'box','on','XScale','log')
xlabel('r','FontSize',font_size)
ylabel('Damon','FontSize',font_size)
title('Error','FontSize',font_size)
axis square

subplot(2,4,6), hold on
plot(diagnostics.D.reg,diagnostics.D.filter_sum,'.','MarkerSize',marker_size2)
plot(diagnostics.D.reg(bf),diagnostics.D.filter_sum(bf),'r.','MarkerSize',marker_size2)
set(gca,'LineWidth',2,'FontSize',20,'box','on','XScale','log')
xlabel('r','FontSize',font_size)
title('\Sigma|K|','FontSize',font_size)
axis square

subplot(2,4,7), hold on
plot(diagnostics.D.reg,diagnostics.D.filter_height,'.','MarkerSize',marker_size2)
plot(diagnostics.D.reg(bf),diagnostics.D.filter_height(bf),'r.','MarkerSize',marker_size2)
set(gca,'LineWidth',2,'FontSize',20,'box','on','XScale','log')
xlabel('r','FontSize',font_size)
title('Filter Height','FontSize',font_size)
axis square

subplot(2,4,8), hold on
plot(diagnostics.D.reg,diagnostics.D.slope,'.','MarkerSize',marker_size2)
plot(diagnostics.D.reg(bf),diagnostics.D.slope(bf),'r.','MarkerSize',marker_size2)
set(gca,'LineWidth',2,'FontSize',20,'box','on','XScale','log')
xlabel('r','FontSize',font_size)
title('Gain','FontSize',font_size)
axis square
