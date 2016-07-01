% compareGainsInFig4
% makes some supplementary plots for fig4.
% this is concerned with various aspects of the variance switching experiment 

function [ax] = compareGainsInFig4(gains,ax,make_plot)

% unpack data
v2struct(gains)

opacity = .5;

if nargin < 2
	figure('outerposition',[0 0 1200 850],'PaperUnits','points','PaperSize',[1200 850]); hold on
	for i = 1:6
		ax(i) = subplot(2,3,i); hold on
	end
	make_plot = true(6,1);
end

% plot transduction gain~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if make_plot(1)
	plot(ax(1),s_hi,gL_hi,'+','Color',[1 opacity opacity])
	errorbar(ax(1),nanmean(s_hi),nanmean(gL_hi),nanstd(gL_hi),'r','LineWidth',4,'Marker','o','MarkerSize',10);

	plot(s_lo,gL_lo,'+','Color',[opacity opacity 1])
	errorbar(nanmean(s_lo),nanmean(gL_lo),nanstd(gL_lo),'b','LineWidth',4,'Marker','o','MarkerSize',10);
	xlabel(ax(1),'\sigma_{Stimulus} (V)')
	ylabel(ax(1),'Transduction gain (mV/V)')
end

% plot firing gain ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if make_plot(2)
	plot(s_hi,gF_hi,'+','Color',[1 opacity opacity])
	errorbar(nanmean(s_hi),nanmean(gF_hi),nanstd(gF_hi),'r','LineWidth',4,'Marker','o','MarkerSize',10);

	plot(s_lo,gF_lo,'+','Color',[opacity opacity 1])
	errorbar(nanmean(s_lo),nanmean(gF_lo),nanstd(gF_lo),'b','LineWidth',4,'Marker','o','MarkerSize',10);
	xlabel(ax(2),'\sigma_{Stimulus} (V)')
	ylabel(ax(2),'Firing gain (Hz/mV)')
end

% plot total gain~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if make_plot(3)
	plot(s_hi,gT_hi,'+','Color',[1 opacity opacity])
	errorbar(nanmean(s_hi),nanmean(gT_hi),nanstd(gT_hi),'r','LineWidth',4,'Marker','o','MarkerSize',10);

	plot(s_lo,gT_lo,'+','Color',[opacity opacity 1])
	errorbar(nanmean(s_lo),nanmean(gT_lo),nanstd(gT_lo),'b','LineWidth',4,'Marker','o','MarkerSize',10);
	xlabel(ax(3),'\sigma_{Stimulus} (V)')
	ylabel(ax(3),'Total gain (Hz/mV)')
end

% plot total gain vs. product of individual gains ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if make_plot(4)
	plot(gL_hi.*gF_hi,gT_hi,'r+')
	plot(gL_lo.*gF_lo,gT_lo,'b+')
	z = max([gL_hi(:); gF_hi(:); gT_hi(:); gT_lo(:); gF_lo(:); gL_lo(:)]);
	plot(ax(4),[0 z],[0 z],'k--')
	xlabel(ax(4),'gain_{LFP} \times gain_{firing} (Hz/V)'); ylabel(ax(4),'Total gain (Hz/V)')
end

% compute fractional gains ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fL = gL_lo./gL_hi;
fF = gF_lo./gF_hi;
fT = gT_lo./gT_hi;

% plot fractional gains vs. product of fractional gains ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if make_plot(5)
	plot(ax(5),fL.*fF,fT,'k+')
	z = max(fF);
	plot(ax(5),[0 z],[0 z],'k--')
	xlabel(ax(5),['Product of fold changes' char(10) 'in gains at firing (\gamma_F) and LFP (\gamma_L)'])
	ylabel(ax(5),['Overall fold change' char(10) 'in gain (\gamma_{Total})'])
end


% compute ratios, contributions, etc ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if make_plot(6)
	clear temp
	errorbar(ax(6),1,mean(fT),std(fT),'LineWidth',4,'Marker','o','MarkerSize',10)
	errorbar(ax(6),2,mean(fL.*fF),std(fL.*fF),'LineWidth',4,'Marker','o','MarkerSize',10)
	errorbar(ax(6),3,mean(fL)*mean(fF),(std(fL)+std(fF))/2,'LineWidth',4,'Marker','o','MarkerSize',10)
	temp(1) = errorbar(4,mean(fL),std(fL),'LineWidth',4,'Marker','o','MarkerSize',10);
	temp(2) = errorbar(5,mean(fF),std(fF),'LineWidth',4,'Marker','o','MarkerSize',10);
	set(ax(6),'XTick',1:6,'XTickLabel',{'\langle\gamma_{Total}\rangle','\langle\gamma_{L} \gamma_{F}\rangle','\langle\gamma_{L}\rangle\langle\gamma_{F}\rangle','\langle\gamma_{L}\rangle','\langle\gamma_{F}\rangle',''},'XTickLabelRotation',45)

	% add a pie chart showing fraction contributed
	ax(7) = axes();
	ax(7).Position(1) = ax(6).Position(1)+ax(6).Position(3)-.03;
	ax(7).Position(2) = ax(6).Position(2)+ax(6).Position(4)/2;
	ax(7).Position(3) = .07;
	ax(7).Position(4) = .15;
	temp2 = mean(log(fL)./log(fL.*fF));
	p = pie([temp2,1-temp2],[0 1]);
	p(1).FaceColor = temp(1).Color;
	p(3).FaceColor = temp(2).Color;

end

% fix some axes
for i = 1:5
	try
		ax(i).YLim(1) = 0;
		ax(i).XLim(1) = 0;
	catch
	end
end


