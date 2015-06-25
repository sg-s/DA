% 
% 
% created by Srinivas Gorur-Shandilya at 10:31 , 21 April 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function [] = MechanismAnalysis_PlotGain(stim,resp,ParadigmNames,paradigm,use_light)

if use_light
	figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
	nplots = 2;
	gainplot = 1;
	phaseplot = 2;
else
	figure('outerposition',[0 0 1400 500],'PaperUnits','points','PaperSize',[1200 800]); hold on
	nplots = 3;
	gainplot = 2;
	phaseplot = 3;
end

c = parula(length(ParadigmNames)+1);
if ~use_light
	subplot(1,nplots,1), hold on
	% show stimulus distributions in all cases. only from 20s-60s 
	
	clear l
	for i = 1:length(ParadigmNames)
		temp = stim(:,paradigm == i);
		temp(1:20e3,:) = []; % throw out the first 20 seconds
		for j = 1:width(temp)
			[y,x] = hist(temp(:,j),50);
			l(i) = plot(x,y,'Color',c(i,:));
		end
	end
	xlabel('Stimulus (V)')
	ylabel('Count')
end


% show the gain changes, if any
ss = 20; % subsample for plot
fmax = 0;
subplot(1,nplots,gainplot), hold on
plot([0 100],[0 100],'k--')

% find which part of the trace we can use
rm_this=(any(isnan(resp')));
resp(rm_this,:) = [];
stim(rm_this,:) = [];

% figure out which is the response with no background
resp0 = resp(:,paradigm == 1);
if width(resp0) > 1
	resp0 = mean2(resp0);
end
L = {};
clear l
rel_gain = ones(length(ParadigmNames),1);
mean_response = NaN(length(ParadigmNames),1);
for i = 2:length(ParadigmNames)
	temp = resp(:,paradigm == i);
	if width(temp)>1
		temp  =mean2(temp);
	end
	plot(resp0(1:ss:end),temp(1:ss:end),'.','Color',c(i,:));
	fmax = max([fmax max(resp0) max(temp)]);

	% fit lines to points
	cf = fit(resp0(:),temp(:),'poly1');
	l(i)=plot(sort(resp0),cf(sort(resp0)),'Color',c(i,:));
	L{i} = strcat('Rel. gain = ',oval(cf.p1));
	rel_gain(i) = cf.p1;
	mean_response(i) = mean(temp);
	mean_response(1) = mean(resp0);

end
set(gca,'XLim',[0 fmax],'YLim',[0 fmax])
l(1) = [];
L(1) = [];
legend(l,L,'Location','southeast')
if use_light
	xlabel('Response to light flicker (Hz)')
	ylabel('Response to light flicker + odour (Hz)')
else
	xlabel('Response to odour flicker (Hz)')
	ylabel('Response to odour flicker + light (Hz)')
end


% make the gain phase plot
alldata.stim = stim;
alldata.resp = resp;
alldata.ParadigmNames = ParadigmNames;
alldata.paradigm = paradigm;
subplot(1,nplots,phaseplot), hold on
if use_light
	GainPhasePlot(alldata,'r')
else
	GainPhasePlot(alldata,'b')
end
