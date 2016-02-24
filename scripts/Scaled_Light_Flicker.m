% 
% 
% created by Srinivas Gorur-Shandilya at 9:30 , 15 February 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

pHeader;

%% Scaled Light Flicker
% In this document, we present scaled Gaussian light inputs to ORNs expressing Chrimson and measuring their firing gain. The following figure shows what the stimulus looks like, and shows the gain as a function of the stimulus mean, for a subset of the data where the contrast is relatively constant. 

p = '/local-data/DA-paper/Chrimson/scaled-light-flicker/v2';
[~, ~, fA, paradigm, orn, fly, AllControlParadigms] = consolidateData(p,true);

x = [.75:0.05:1.1 1.2:.1:3.6 3.8 4 4.2 4.5];
x = [0:0.1:0.7 x];
y = [0 4 12.8 24.1 36.2 48.1 60 70 95 114 134 151 167 184 201 219 236 252 269 283 299 314 328 343 357 370 383 394 404 413 419 422 422 421 419 418 417];
y = [0*(0:0.1:0.7) y ];

light_power_fit = fit(x(:),y(:),'spline');

% make new vectors for mean and range of the stimulus
a = 25e3; z = 55e3;
LED = NaN*fA;
LED_voltage = LED;
light_s = NaN*paradigm; light_m = NaN*paradigm;
for i = 1:length(paradigm)
	temp = AllControlParadigms(paradigm(i)).Name;
	LED_voltage(:,i) = (AllControlParadigms(paradigm(i)).Outputs(1,1:10:end));
	LED(:,i) = light_power_fit(AllControlParadigms(paradigm(i)).Outputs(1,1:10:end));
end

% remove some junk
rm_this = sum(fA) == 0;
LED(:,rm_this) = [];
fA(:,rm_this) = [];
paradigm(rm_this) = [];
orn(rm_this) = [];
fly(rm_this) = [];
LED_voltage(:,rm_this) = [];

% extract filters
[K,fp,gain,gain_err] = extractFilters(LED,fA,'use_cache',true,'a',a,'z',z);


% show the stimulus distributions
figure('outerposition',[0 0 1200 800],'PaperUnits','points','PaperSize',[1300 800]); hold on
for i = 1:6
	ax(i) = subplot(2,3,i); hold on
end
light_levels = unique(mean(LED(a:z,:)));
all_light_levels = (mean(LED(a:z,:)));
c = parula(length(light_levels)+1);
hx = logspace(0,log10(500),200);
for i = 4:length(c)-1
	these_paradigms = find(all_light_levels == light_levels(i));

	% plot the mean vs. the std dev of the light stimulus
	plot(ax(1),mean(mean(LED(a:z,these_paradigms))),mean(std(LED(a:z,these_paradigms))),'+','Color',c(i,:))

	% also plot the contrast
	plot(ax(2),mean(mean(LED(a:z,these_paradigms))),mean(std(LED(a:z,these_paradigms)))./mean(mean(LED(a:z,these_paradigms))),'+','Color',c(i,:))

	temp = LED(a:z,these_paradigms(1));
	temp(temp<0) = [];
	y = histcounts(temp,hx);
	y = y/sum(y);
	xx = hx(1:end-1) + mean(diff(hx));
	plot(ax(3),xx,y,'Color',c(i,:))

	% plot the filters
	tK = nanmean(K(:,paradigm==i),2);
	filtertime = -99:600;
	plot(ax(4),filtertime,tK,'Color',c(i,:))

	x = fp(:,these_paradigms);
	if width(x) > 1
		x = nanmean(x,2);
	end
	y = fA(:,these_paradigms);
	if width(y) > 1
		y = nanmean(y,2);
	end
	axes(ax(5))
	plotPieceWiseLinear(x,y,'Color',c(i,:),'nbins',40);
	plot(ax(6),all_light_levels(these_paradigms),gain(these_paradigms),'+','Color',c(i,:))
end

% fo = fitoptions('power1');
% fo.Lower = [-Inf -1]; fo.Upper = [Inf -1];
these_paradigms = (all_light_levels > light_levels(3));

ff = fit(all_light_levels(these_paradigms)',gain(these_paradigms),'power1');
fx = [min(all_light_levels(these_paradigms)) max(all_light_levels(these_paradigms))];
l = plot(ax(6),fx,ff(fx),'k--');
legend(l,['\alpha = ' oval(ff.b)])

xlabel(ax(1),'\mu_{Light} (\muW)')
ylabel(ax(1),'\sigma_{Light} (\muW)')

xlabel(ax(2),'\mu_{Light} (\muW)')
ylabel(ax(2),'\sigma/\mu')
set(ax(2),'YLim',[0 .5])

xlabel(ax(3),'Light Power (\muW)')
ylabel(ax(3),'Probability')
set(ax(3),'XScale','log')

xlabel(ax(4),'Filter Lag (ms)')
ylabel(ax(4),'Filter ')

set(ax(5),'XScale','log')
xlabel(ax(5),'Projected Stimulus (\muW)')
ylabel(ax(5),'ORN Response (Hz)')
xlabel(ax(6),'Mean Stimulus (\muW)')
set(ax(6),'YScale','log','XScale','log')
ylabel(ax(6),'Gain (Hz/\muW)')

prettyFig('FixLogY=true;','FixLogX=true;','fs=18;')

if being_published	
	snapnow	
	delete(gcf)
end

%%
% We now show the full dataset, including some parts of the data where the stimulus is very low, where, due to some yet-to-be-determined error, the contrast drops.  

figure('outerposition',[0 0 1200 800],'PaperUnits','points','PaperSize',[1300 800]); hold on
for i = 1:6
	ax(i) = subplot(2,3,i); hold on
end
light_levels = unique(mean(LED(a:z,:)));
all_light_levels = (mean(LED(a:z,:)));
c = parula(length(light_levels)+1);
hx = logspace(0,log10(500),200);
for i = 1:length(c)-1
	these_paradigms = find(all_light_levels == light_levels(i));

	% plot the mean vs. the std dev of the light stimulus
	plot(ax(1),mean(mean(LED(a:z,these_paradigms))),mean(std(LED(a:z,these_paradigms))),'+','Color',c(i,:))

	% also plot the contrast
	plot(ax(2),mean(mean(LED(a:z,these_paradigms))),mean(std(LED(a:z,these_paradigms)))./mean(mean(LED(a:z,these_paradigms))),'+','Color',c(i,:))

	temp = LED(a:z,these_paradigms(1));
	temp(temp<0) = [];
	y = histcounts(temp,hx);
	y = y/sum(y);
	xx = hx(1:end-1) + mean(diff(hx));
	plot(ax(3),xx,y,'Color',c(i,:))

	% plot the filters
	tK = nanmean(K(:,paradigm==i),2);
	filtertime = -99:600;
	plot(ax(4),filtertime,tK,'Color',c(i,:))

	x = fp(:,these_paradigms);
	if width(x) > 1
		x = nanmean(x,2);
	end
	y = fA(:,these_paradigms);
	if width(y) > 1
		y = nanmean(y,2);
	end
	axes(ax(5))
	plotPieceWiseLinear(x,y,'Color',c(i,:),'nbins',40);
	plot(ax(6),all_light_levels(these_paradigms),gain(these_paradigms),'+','Color',c(i,:))
end


these_paradigms = (all_light_levels >= light_levels(1));

ff = fit(all_light_levels(these_paradigms)',gain(these_paradigms),'power1');
fx = [min(all_light_levels(these_paradigms)) max(all_light_levels(these_paradigms))];
l = plot(ax(6),fx,ff(fx),'k--');
legend(l,['\alpha = ' oval(ff.b)])

xlabel(ax(1),'\mu_{Light} (\muW)')
ylabel(ax(1),'\sigma_{Light} (\muW)')

xlabel(ax(2),'\mu_{Light} (\muW)')
ylabel(ax(2),'\sigma/\mu')
set(ax(2),'YLim',[0 .5])

xlabel(ax(3),'Light Power (\muW)')
ylabel(ax(3),'Probability')
set(ax(3),'XScale','log')

xlabel(ax(4),'Filter Lag (ms)')
ylabel(ax(4),'Filter ')

set(ax(5),'XScale','log')
xlabel(ax(5),'Projected Stimulus (\muW)')
ylabel(ax(5),'ORN Response (Hz)')
xlabel(ax(6),'Mean Stimulus (\muW)')
set(ax(6),'YScale','log','XScale','log')
ylabel(ax(6),'Gain (Hz/\muW)')

prettyFig('FixLogY=true;','FixLogX=true;','fs=18;')

if being_published	
	snapnow	
	delete(gcf)
end


%% Kinetics of Response as a function of mean stimulus
% In this section, we determine if the firing  responses slow down or speed up as a function of the mean stimulus. We do so buy computing the cross-correlation functions between the stimulus and the response for each trial.

firing_xcorr = NaN(1e3,width(LED));
firing_peaks = NaN*paradigm;
firing_xcorr_max = NaN*paradigm;

for i = 1:width(LED)
	s = LED(25e3:45e3,i)-mean(LED(25e3:45e3,i)); s = s/std(s);
	f = fA(25e3:45e3,i)-mean(fA(25e3:45e3,i)); f = f/std(f);

	temp = xcorr(f,s);
	firing_xcorr(:,i) = temp(19.5e3+1:20.5e3);
	[firing_xcorr_max(i),firing_peaks(i)] = max(firing_xcorr(:,i));
	firing_xcorr(:,i) = firing_xcorr(:,i)/max(firing_xcorr(:,i));
end
firing_peaks = firing_peaks - 500;
firing_peaks(firing_peaks<0) = NaN;

figure('outerposition',[0 0 1200 400],'PaperUnits','points','PaperSize',[1300 400]); hold on
ax(1) = subplot(1,3,1); hold on
plot([0 0],[-1 1],'k--')
title('Firing xcorr functions')
ylabel('Cross correlation (norm)')
xlabel('Lag (ms)')
for i = 1:length(c)- 1
	try
		plot(ax(1),-499:500,firing_xcorr(:,paradigm==i),'Color',c(i,:),'LineWidth',3)
	catch
	end
end
set(ax(1),'XLim',[-10 250])
subplot(1,3,2), hold on
plot(abs(firing_xcorr_max)/20e3,firing_peaks,'k+')
xlabel('Absolute Peak correlation')
ylabel('Location of peak (ms)')

subplot(1,3,3), hold on
l = plot(nanmean(LED(35e3:45e3,:)),firing_peaks,'k+');
xlabel('Mean Light Power (\muW)')
ylabel('Delay (ms)')
legend(l,['\rho = ' oval(spear(nanmean(LED(35e3:45e3,:))',firing_peaks(:)))])


prettyFig('plw=1.3;','lw=1.5;','fs=14;','FixLogX=true;','FixLogY=0;')

if being_published
	snapnow
	delete(gcf)
end

%% Version Info
%
pFooter;


