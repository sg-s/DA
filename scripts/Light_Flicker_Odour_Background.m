% Light_Flicker_Odour_Background.m
% 
% created by Srinivas Gorur-Shandilya at 4:06 , 08 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


pHeader;

%% Light Flicker Odour Background
% This document looks at data where we present a odour background and a light flicker, and looks at how (and if) the gain to the light changes with odour background.

x = [.75:0.05:1.1 1.2:.1:3.6 3.8 4 4.2 4.5];
x = [0:0.1:0.7 x];
y = [0 4 12.8 24.1 36.2 48.1 60 70 95 114 134 151 167 184 201 219 236 252 269 283 299 314 328 343 357 370 383 394 404 413 419 422 422 421 419 418 417];
y = [0*(0:0.1:0.7) y ];

led_power_func = fit(x(:),y(:),'spline');

if ~exist('od','var')
	p = '/local-data/DA-paper/Chrimson/light-odor-combinations/light-flicker/v3';
	od = raw2ORNData(p,'led_power_func',led_power_func,'use_led',true,'filter_LFP',false);

	uts = false(length(od(1,1).stimulus),1);
	uts(10e3:end-10e3) = true;
	for i = 1:size(od,1)
		for j = 1:size(od,2)
			if od(i,j).n_trials > 0
				disp([i j])
				od(i,j).use_this_segment = uts;
				od(i,j) = backOutFilters(od(i,j));
			end
		end
	end

	% also create another object with the odour values
	od_odour = raw2ORNData(p,'led_power_func',led_power_func,'use_led',false,'filter_LFP',false);

end

% figure out paradigm ordering 
mean_s = NaN(size(od_odour,2),1);
for i = 1:size(od_odour,2)
	mean_s(i) = nanmean(nanmean([od_odour(:,i).stimulus]));
end
[~,idx] = sort(mean_s);
for i = 1:size(od_odour,1)
	od_odour(i,:) = permute(od_odour(i,:),idx);
	od(i,:) = permute(od(i,:),idx);
end

% remove the baseline PID from all the PID values
PID_baseline = nanmean(nanmean([od_odour(:,1).stimulus]));
for i = 1:size(od_odour,1)
	od_odour(i,1).stimulus = 0*od_odour(i,1).stimulus + 1e-4; % we will reset this to 0 in the plot
	for j = 2:size(od_odour,2)
		od_odour(i,j).stimulus = od_odour(i,j).stimulus - PID_baseline;
	end
end

%%
% First, we back out filters and plot the nonlinearities in the light flicker regime, and plot them colour coded by odour background level.

figure('outerposition',[0 0 1600 800],'PaperUnits','points','PaperSize',[1600 800]); hold on
clear ax
for i = size(od,1)*3:-1:1
	ax(i) = subplot(3,size(od,1),i); hold on
end
filtertime = 1e-3*(1:700) - .1;
for i = 1:size(od,1)
	c = parula(length(find([od(i,:).n_trials]))+1);
	ci = 1;
	gain = NaN(size(od,2),1);
	odour_levels = NaN(size(od,2),1);
	for j = 1:size(od,2)
		if od(i,j).n_trials > 0
			K = (nanmean(od(i,j).K_firing,2));
			x = (nanmean(od(i,j).firing_projected(uts,:),2)); x = x(:);
			y = (nanmean(od(i,j).firing_rate(uts,:),2)); y = y(:);
			plot(ax(i),filtertime,K,'Color',c(ci,:));
			axes(ax(i+size(od,1)))
			[~,data] = plotPieceWiseLinear(x,y,'nbins',50,'Color',c(ci,:));
			middle_segment = (y > max(y)/3 & y < 2*max(y)/3);
			ff = fit(x(middle_segment),y(middle_segment),'poly1');
			gain(j) = ff.p1;
			odour_levels(j) = nanmean((nanmean(od_odour(i,j).stimulus(:,:),2)));
			ci = ci+1;
		end
	end
	plot(ax(i+size(od,1)*2),odour_levels,gain,'k+')
end

% equalise axes everywhere
set(ax(size(od,1)*2+1:end),'YLim',[0 1.5])

% add labels, etc
for i = 1:size(od,1)
	ylabel(ax(i),'Filter')
	xlabel(ax(i),'Filter lag (s)')
	title(ax(i),['ORN #' oval(i)])
end

for i = size(od,1)+1:size(od,1)*2
	ylabel(ax(i),'ORN Response (Hz)')
	xlabel(ax(i),'Projected Stimulus (\muW)')
end

for i = size(od,1)*2+1:length(ax)
	ylabel(ax(i),'Gain (Hz/\muW)')
	xlabel(ax(i),'Odor Stimulus (V)')
end

prettyFig('fs=18;')

if being_published	
	snapnow	
	delete(gcf)
end

%%
% Now, we summarise how the gain changes in all the neurons in one figure:

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
for i = 1:size(od,1)
	c = lines(size(od,1));
	clear gain odour_levels
	gain = NaN(size(od,2),1);
	odour_levels = NaN(size(od,2),1);
	for j = 1:size(od,2)
		if od(i,j).n_trials > 0
			K = (nanmean(od(i,j).K_firing,2));
			x = (nanmean(od(i,j).firing_projected(uts,:),2)); x = x(:);
			y = (nanmean(od(i,j).firing_rate(uts,:),2)); y = y(:);

			[~,data] = plotPieceWiseLinear(x,y,'nbins',50,'make_plot',false);
			middle_segment = (y > max(y)/3 & y < 2*max(y)/3);
			ff = fit(x(middle_segment),y(middle_segment),'poly1');
			gain(j) = ff.p1;
			odour_levels(j) = nanmean((nanmean(od_odour(i,j).stimulus(:,:),2)));
		end
	end
	rm_this = isnan(odour_levels) | isnan(gain);
	odour_levels(rm_this) = []; gain(rm_this) = [];
	plot(odour_levels,gain,'-+','Color',c(i,:))
end

set(gca,'XScale','log','YLim',[0 1])
xlabel('Odour Stimulus (V)')
ylabel('Gain (Hz/\muW)')
prettyFig('fs=18;')

if being_published	
	snapnow	
	delete(gcf)
end




%%
% Now, we show that the background odour activates the ORN. In our experiments, the odour is turned on 5 seconds before the light starts to fluctuate. We plot the ORN response during this initial period as a function of odour concentration. 

figure('outerposition',[100 100 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
c = lines(size(od,1));
for i = 1:size(od,1)
	% for each neuron
	s = NaN(size(od,2),1);
	f = NaN(size(od,2),1);
	for j = 1:size(od,2)
		% for each paradigm
		f(j) = nanmean(nanmean(od(i,j).firing_rate(1:5e3,:)));
		s(j) = nanmean(nanmean(od_odour(i,j).stimulus(1:5e3,:)));
	end
	rm_this = isnan(s) | isnan(f);
	s(rm_this) = []; f(rm_this) = [];
	plot(s,f,'-+','Color',c(i,:))
end
set(gca,'XScale','log')
xlabel('Mean Odour Stimulus (V)')
ylabel('Odour-induced Firing rate (Hz)')

prettyFig('FixLogX=true;')

if being_published	
	snapnow	
	delete(gcf)
end





%% Version Info
%
pFooter;


