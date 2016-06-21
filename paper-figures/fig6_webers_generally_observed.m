
%% Weber's Law generally observed
% In this document we look at a number of odorants and receptors and check if Weber's Law is true in the LFP and the firing rate. We do so by presenting Gaussian stimuli with increasing means. 


pHeader;
dm = dataManager;

% define what we want to work on
data_hashes = {'bcd4cf4fe12817d084a2b06f981161ee','cd6753c0e4cf02895cd5e2c5cb58aa1a','3ea08ccfa892c6545d74bbdaaa6cbee1','a33723c87a1216b274750734b4ee8820'};
odour_names = {'1-pentanol','1-pentanol','2-butanone','isoamyl-acetate'};
orn_names = {'ab3A','ab2A','ab2A','pb1A'};

f1 = figure('outerposition',[0 0 1400 700],'PaperUnits','points','PaperSize',[1400 700]); hold on
clear ax axs
for i = 1:2*length(data_hashes)
	ax(i) = subplot(2,length(data_hashes),i); hold on
end

f2 = figure('outerposition',[0 0 1400 700],'PaperUnits','points','PaperSize',[1400 700]); hold on
for i = 1:2*length(data_hashes)
	axs(i) = subplot(2,length(data_hashes),i); hold on
end

% core loop
for i = length(data_hashes):-1:1
	clear cdata
	cdata = consolidateData2(dm.getPath(data_hashes{i}));
	cdata = cleanMSGdata(cdata);

	% plot gain as we normally calculate it
	plot_here = [ax(i) ax(i+length(data_hashes))];
	plotMSGGain(cdata,plot_here);

	% for the supp. figure, plot gain as ratios of variances
	plot_here = [axs(i) axs(i+length(data_hashes))];
	plotMSGGain2(cdata,plot_here);
end

% add titles 
for i = 1:length(data_hashes)
	t = [orn_names{i} char(10) odour_names{i}];
	title(ax(i),t);
	title(axs(i),t);
end

prettyFig(f1,'FixLogX',true);
prettyFig(f2,'FixLogX',true);

for i = 1:length(ax)
	deintersectAxes(ax(i));
	deintersectAxes(axs(i));
end

if being_published
	snapnow
	delete(gcf)
end




%% Version Info
%
pFooter;
