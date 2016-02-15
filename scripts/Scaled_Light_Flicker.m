% 
% 
% created by Srinivas Gorur-Shandilya at 9:30 , 15 February 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

pHeader;

%% Scaled Light Flicker
% In this document, we present scaled Gaussian light inputs to ORNs expressing Chrimson and measuring their firing gain. The following figure shows what the stimulus looks like, and shows the gain as a function of the stimulus mean. 

p = '/local-data/DA-paper/Chrimson/scaled-light-flicker';
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

% because we fucked up in constructing the stimuus, only keep the paradigms where the fuck up isn't so bad
ok_paradigms = ((mean((LED_voltage(a:z,:))<.8)) < .1);
rm_this = ~ok_paradigms | sum(fA) == 0;
LED(:,rm_this) = [];
LED_voltage(:,rm_this) = [];
fA(:,rm_this) = [];
paradigm(rm_this) = [];
fly(rm_this) = [];
orn(rm_this) = [];

% extract filters
[K,fp,gain,gain_err] = extractFilters(LED,fA,'use_cache',true,'a',a,'z',z);

% show the stimulus distributions
figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
ax(1) = subplot(1,3,1); hold on
ax(2) = subplot(1,3,2); hold on
ax(3) = subplot(1,3,3); hold on
light_levels = unique(mean(LED(a:z,:)));
all_light_levels = (mean(LED(a:z,:)));
c = parula(length(light_levels)+1);
for i = 1:length(c)-1
	these_paradigms = find(all_light_levels == light_levels(i));
	temp = LED(a:z,these_paradigms(1));
	temp(temp<10) = [];
	[y,x] = histcounts(temp,50);
	y = y/sum(y);
	x = x(1:end-1) + mean(diff(x));
	plot(ax(1),x,y,'Color',c(i,:))

	x = fp(:,these_paradigms);
	if width(x) > 1
		x = nanmean(x,2);
	end
	y = fA(:,these_paradigms);
	if width(y) > 1
		y = nanmean(y,2);
	end
	axes(ax(2))
	plotPieceWiseLinear(x,y,'Color',c(i,:),'nbins',40);
	plot(ax(3),all_light_levels(these_paradigms),gain(these_paradigms),'+','Color',c(i,:))
end
xlabel(ax(2),'Projected Stimulus (\muW)')
ylabel(ax(2),'ORN Response (Hz)')
xlabel(ax(3),'Mean Stimulus (\muW)')
ylabel(ax(3),'Gain (Hz/\muW)')
xlabel(ax(1),'Light Power (\muW)')
ylabel(ax(1),'Probability')

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end


%% Version Info
%
pFooter;


