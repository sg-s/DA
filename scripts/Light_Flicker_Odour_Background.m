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
			od(i,j).use_this_segment = uts;
			od(i,j) = backOutFilters(od(i,j));
		end
	end

	% also create another object with the odour values
	od_odour = raw2ORNData(p,'led_power_func',led_power_func,'use_led',false,'filter_LFP',false);
end

%%
% First, we show that the background odour activates the ORN. In our experiments, the odour is turned on 5 seconds before the light starts to fluctuate. We plot the ORN response during this initial period as a function of odour concentration. 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ax(1) = subplot(1,2,1); hold on
ax(2) = subplot(1,2,2); hold on
for i = size(od,1)
	% for each neuron
	for j = 1:size(od,2)
		% for each paradigm
		f = nanmean(od(i,j).firing_rate(1:5e3,:));
		x = od(i,j).LFP(1,:) - nanmean(od(i,j).LFP(1:5e3,:));
		s = nanmean(od_odour(i,j).stimulus(1:5e3,:));
		plot(ax(1),s,x,'k+')
		plot(ax(2),s,f,'k+')
	end
end
xlabel(ax(1),'Mean Odour Stimulus (V)')
xlabel(ax(2),'Mean Odour Stimulus (V)')

ylabel(ax(1),'Mean \DeltaLFP (mV)')
ylabel(ax(2),'Mean Firing rate (Hz)')

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end


% figure out paradigm ordering 
mean_s = NaN(size(od_odour,2),1);
for i = 1:size(od_odour,2)
	mean_s(i) = nanmean(nanmean(od_odour(:,i).stimulus));
end
[~,idx] = sort(mean_s);
for i = 1:size(od_odour,1)
	od_odour(i,:) = permute(od_odour(i,:),idx);
	od(i,:) = permute(od(i,:),idx);
end

%%
% Now, we back out filters and plot the nonlinearities in the light flicker regime, and plot them colour coded by odour background level.

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ax(1) = subplot(1,2,1); hold on
ax(2) = subplot(1,2,2); hold on
c = parula(size(od,2)+1);
filtertime = 1e-3*(1:700) - .1;
for i = 1:size(od,2)
	K = (nanmean([od(:,i).K_firing],2));
	x = (nanmean([od(:,i).firing_projected(uts,:)],2));
	y = (nanmean([od(:,i).firing_rate(uts,:)],2));
	plot(ax(1),filtertime,K,'Color',c(i,:));
	axes(ax(2))
	plotPieceWiseLinear(x,y,'nbins',50,'Color',c(i,:));
end
xlabel(ax(1),'Fitler Lag (s)')
ylabel(ax(1),'Filter')
xlabel(ax(2),'Projected Stimulus (\muW)')
ylabel(ax(2),'Firing rate (Hz) ')

prettyFig('fs=18;')

if being_published	
	snapnow	
	delete(gcf)
end

%% Version Info
%
pFooter;


