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
	mplots = 1;
	filterplot = 1;
	gainplot = 2;
else
	figure('outerposition',[0 0 1200 800],'PaperUnits','points','PaperSize',[1500 500]); hold on
	nplots = 3;
	mplots = 2;
	filterplot = 2;
	gainplot = 3;
end

c = parula(length(ParadigmNames)+1);
if ~use_light
	subplot(mplots,nplots,1), hold on
	% show stimulus distributions in all cases. 
	
	clear l
	for i = 1:length(ParadigmNames)
		temp = stim(:,paradigm == i);
		temp(1:1e4,:) = [];
		for j = 1:width(temp)
			[y,x] = hist(temp(:,j),50);
			l(i) = plot(x,y,'Color',c(i,:));
		end
	end
	xlabel('Stimulus (V)')
	ylabel('Count')
end

% show filters for each case
subplot(mplots,nplots,filterplot), hold on
clear l
for i = 1:width(stim)
	this_stim = stim(:,i);
	this_resp = resp(:,i);
	this_stim(1:1e4,:) = [];
	this_resp(1:1e4,:) = [];
	[this_K, ~, filtertime_full] = FindBestFilter(this_stim,this_resp,[],'regmax=1;','regmin=1;','filter_length=1999;','offset=500;');
	filtertime_full = filtertime_full*1e-3;
	filtertime = 1e-3*(-200:900);
	this_K = interp1(filtertime_full,this_K,filtertime);
	this_K = this_K/max(this_K);
	l(paradigm(i)) = plot(filtertime,this_K,'Color',c(paradigm(i),:));
end
legend(l,ParadigmNames,'Location','southoutside')
xlabel('Filter Lag (s)')
ylabel('Filter Amplitude (norm)')

% show the gain changes, if any
ss = 20;
fmax = 0;
subplot(mplots,nplots,gainplot), hold on
plot([0 100],[0 100],'k--')
resp0 = resp(:,paradigm == 1);
resp0(1:1e4,:) = [];
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
	temp(1:1e4,:) = [];
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

% compute correlation functions
subplot(mplots,nplots,4), hold on
resp0 = resp(:,paradigm == 1);
if width(resp0) > 1
	resp0 = mean2(resp0);
end
resp0 = resp0 - mean(resp0);
resp0 = resp0/std(resp0);
peak_xcorr = NaN(length(ParadigmNames)-1,1);
for i = 1:length(peak_xcorr)
	this_resp = resp(:,paradigm == (1 + i));
	if width(this_resp) > 1
		this_resp = mean2(this_resp);
	end
	this_resp = this_resp - mean(this_resp);
	this_resp = this_resp/std(this_resp);
	peak_xcorr(i) = finddelay(resp0,this_resp);
	x = xcorr(this_resp(1e4:end),resp0(1e4:end));
	x = x/max(x);
	tx = 1e-3*(1:length(x));
	tx = tx - mean(tx);
	plot(tx,x,'Color',c(i+1,:))
end
set(gca,'XLim',[-.2 .5])
ylabel('Cross Correlation')
xlabel('Lag (s)')


% plot peak vs gain
subplot(mplots,nplots,5), hold on
plot(rel_gain(2:end),peak_xcorr,'k+')
xlabel('Relative Gain')
ylabel('Delay (ms)')

% plot gain vs response
subplot(mplots,nplots,6), hold on
plot(mean_response,rel_gain,'k+')
ylabel('Relative Gain')
xlabel('Mean Response (Hz)')
