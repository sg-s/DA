% fig_light_mechanism.m
% 
% created by Srinivas Gorur-Shandilya at 9:50 , 28 December 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


%% Mechanism of Gain Control: Optogenetic Activation

pHeader;

% build the LED map
x = [.75:0.05:1.1 1.2:.1:3.6 3.8 4 4.2 4.5];
x = [0:0.1:0.7 x];
y = [0 4 12.8 24.1 36.2 48.1 60 70 95 114 134 151 167 184 201 219 236 252 269 283 299 314 328 343 357 370 383 394 404 413 419 422 422 421 419 418 417];
y = [0*(0:0.1:0.7) y ];
light_power_fit = fit(x(:),y(:),'spline');

figure('outerposition',[0 0 850 850],'PaperUnits','points','PaperSize',[850 850]); hold on

% ##       ####  ######   ##     ## ######## 
% ##        ##  ##    ##  ##     ##    ##    
% ##        ##  ##        ##     ##    ##    
% ##        ##  ##   #### #########    ##    
% ##        ##  ##    ##  ##     ##    ##    
% ##        ##  ##    ##  ##     ##    ##    
% ######## ####  ######   ##     ##    ##    

%  ######   #######  ##    ## ######## ########     ###     ######  ######## 
% ##    ## ##     ## ###   ##    ##    ##     ##   ## ##   ##    ##    ##    
% ##       ##     ## ####  ##    ##    ##     ##  ##   ##  ##          ##    
% ##       ##     ## ## ## ##    ##    ########  ##     ##  ######     ##    
% ##       ##     ## ##  ####    ##    ##   ##   #########       ##    ##    
% ##    ## ##     ## ##   ###    ##    ##    ##  ##     ## ##    ##    ##    
%  ######   #######  ##    ##    ##    ##     ## ##     ##  ######     ##    


p = '/local-data/DA-paper/Chrimson/MSG/ok';
[~, ~, fA, paradigm, orn, fly, AllControlParadigms] = consolidateData(p,true);


% make new vectors for mean and range of the stimulus
a = 25e3; z = 55e3;
LED = NaN*fA;
r = NaN*paradigm; m = NaN*paradigm;
light_s = NaN*paradigm; light_m = NaN*paradigm;
for i = 1:length(paradigm)
	temp = AllControlParadigms(paradigm(i)).Name;
	r(i) = (str2double(temp(3:strfind(temp,'_')-1)));
	m(i) = (str2double(temp(strfind(temp,'_')+3:end)));
	LED(:,i) = light_power_fit(AllControlParadigms(paradigm(i)).Outputs(1,1:10:end));
	light_s(i) = std(LED(a:z,i));
	light_m(i) = mean(LED(a:z,i));
end

% extract filters
[K,fp,gain,gain_err] = extractFilters(LED,fA,'use_cache',true,'a',a,'z',z);

contrast_levels = sort(unique(light_s(m==2)));
light_m = NaN*contrast_levels;
c = parula(length(contrast_levels)+1);
for i = 1:length(contrast_levels)
	light_m(i) = mean(mean(LED(a:z,light_s == contrast_levels(i))));
	y = (fA(a:z,light_s == contrast_levels(i)));
	if ~isvector(y)
		y = nanmean(y,2);
	end
	x = (fp(a:z,light_s == contrast_levels(i)));
	if ~isvector(x)
		x = nanmean(x,2);
	end
	x = x - mean(x) + light_m(i);
	
	[yy,xx] = histcounts(x,0:600);
	yy = yy/sum(yy);
end

subplot(3,3,8), hold on
% plot i/o curves and distributions for lowest contrast
do_these = [4 1];
c = [0 0 1; 0 0 0; 0 0 0 ; 1 0 0];
for i = do_these
	y = (fA(a:z,light_s == contrast_levels(i)));
	if ~isvector(y)
		y = nanmean(y,2);
	end
	x = (fp(a:z,light_s == contrast_levels(i)));
	if ~isvector(x)
		x = nanmean(x,2);
	end
	x = x - mean(x) + mean(mean(LED(a:z,light_s == contrast_levels(i))));

	[~,data] = plotPieceWiseLinear(x,y,'make_plot',false,'nbins',40);
	plot(data.x,data.y,'Color',c(i,:),'LineWidth',3)
end
xlabel('Projected Stimulus (\muW)')
ylabel('Firing Rate (Hz)')

% plot gain vs. contrast
subplot(3,3,9), hold on
plot(light_s(m==2),gain(m==2),'k+')
xlabel('\sigma_{Light} (\muW)')
ylabel('ORN Gain (Hz/\muW)')


%    ##       ####  ######   ##     ## ########    ########   ######   
%    ##        ##  ##    ##  ##     ##    ##       ##     ## ##    ##  
%    ##        ##  ##        ##     ##    ##       ##     ## ##        
%    ##        ##  ##   #### #########    ##       ########  ##   #### 
%    ##        ##  ##    ##  ##     ##    ##       ##     ## ##    ##  
%    ##        ##  ##    ##  ##     ##    ##       ##     ## ##    ##  
%    ######## ####  ######   ##     ##    ##       ########   ######   


p = '/local-data/DA-paper/Chrimson/light-odor-combinations/odour-flicker/';
[PID, LFP, fA, paradigm, orn, fly, AllControlParadigms, paradigm_hashes] = consolidateData(p,true);
a = 15e3; z = 55e3;

% clean up the data
fA(:,sum(fA) == 0) = NaN;

% figure out which paradigms to plot
light_background_paradigms = false(length(paradigm),1);
odour_background_paradigms = false(length(paradigm),1);
for i = 1:length(paradigm)
	if any(strfind(AllControlParadigms(paradigm(i)).Name,'0%')) || any(strfind(AllControlParadigms(paradigm(i)).Name,'Light'))
		light_background_paradigms(i) = true;
	end
	if any(strfind(AllControlParadigms(paradigm(i)).Name,'%')) 
		odour_background_paradigms(i) = true;
	end
end

% calculate the mean light power in every trial
mean_light_power = zeros(length(paradigm),1);
mean_PID = mean(PID(a:z,:));
for i = 1:length(paradigm)
	this_p = paradigm(i);
	try
		mean_light_power(i) = light_power_fit(str2double(AllControlParadigms(this_p).Name(7:strfind(AllControlParadigms(this_p).Name,'V')-1)));
	catch
	end
end
mean_light_power(isnan(mean_light_power)) = 0;

% filter all LFPs
filtered_LFP = LFP;
% filter the LFP
for i = 1:width(LFP)
	filtered_LFP(:,i) = LFP(:,i) - fastFiltFilt(ones(1e4,1),1e4,LFP(:,i));
	filtered_LFP(:,i) = filtered_LFP(:,i)*10; % to get the units right, now in mV
end


% extract filters everywhere for LFP and gain
[K_LFP,LFP_pred,LFP_gain] = extractFilters(PID,filtered_LFP,'use_cache',true,'a',a,'z',z);
[K_fA,fA_pred,fA_gain] = extractFilters(PID,fA,'use_cache',true,'a',a,'z',z);


% % plot the firing gain for odour and odour
% plot_these_paradigms = sort(unique(paradigm(odour_background_paradigms)));
% x = []; y = []; base_gain = [];
% for i = 1:length(paradigm)
% 	if odour_background_paradigms(i)
% 		ci = find(paradigm(i) == plot_these_paradigms);
% 		plot(axes_handles(7),mean_PID(i),fA_gain(i),'k+');
% 		x = [x mean_PID(i)];
% 		y = [y fA_gain(i)];
% 		if ci == 1
% 			base_gain = [base_gain fA_gain(i)];
% 		end
% 	end
% end
% set(axes_handles(7),'XScale','log','YScale','log')
% rm_this = isnan(x) | isnan(y);
% x(rm_this) = []; y(rm_this) = [];
% ff = fit(x(:),y(:),'power1','Lower',[-Inf -1],'Upper',[Inf -1]);
% % also plot a flat line for comparison
% plot(axes_handles(7),[min(x) max(x)],[nanmean(base_gain) nanmean(base_gain)],'k--')
% plot(axes_handles(7),sort(x),ff(sort(x)),'k')
% ylabel(axes_handles(7),'ORN Gain (Hz/V)')
% xlabel(axes_handles(7),'Odor concentration (V)')


% % also plot the LFP gain
% x = []; y = []; base_gain = [];
% for i = 1:length(paradigm)
% 	if odour_background_paradigms(i)
% 		ci = find(paradigm(i) == plot_these_paradigms);
% 		plot(axes_handles(8),mean_PID(i),LFP_gain(i),'k+');
% 		x = [x mean_PID(i)];
% 		y = [y LFP_gain(i)];
% 		if ci == 1
% 			base_gain = [base_gain LFP_gain(i)];
% 		end
% 	end
% end
% plot(axes_handles(8),[min(x) max(x)],[nanmean(base_gain) nanmean(base_gain)],'k--')
% set(axes_handles(8),'XScale','log','YScale','log','YLim',[.05 2])
% rm_this = isnan(x) | isnan(y);
% x(rm_this) = []; y(rm_this) = [];
% ff = fit(x(:),y(:),'power1','Lower',[-Inf -1],'Upper',[Inf -1]);
% plot(axes_handles(8),sort(x),ff(sort(x)),'k')
% ylabel(axes_handles(8),'LFP Gain (mV/V)')
% xlabel(axes_handles(8),'Odor concentration (V)')

% plot the light background odour flicker cases -- firing gain
plot_these_paradigms = sort(unique(paradigm(light_background_paradigms)));
subplot(3,3,3), hold on
x = []; y = []; base_gain = [];
for i = 1:length(paradigm)
	if light_background_paradigms(i)
		ci = find(paradigm(i) == plot_these_paradigms);
		plot(mean_light_power(i),fA_gain(i),'k+');
		x = [x; mean_light_power(i)];
		y = [y; fA_gain(i)];
		if ci == 1
			base_gain = [base_gain fA_gain(i)];
		end
	end
end

plot([min(x) max(x)],[nanmean(base_gain) nanmean(base_gain)],'k--')
set(gca,'XScale','linear','YScale','log','YLim',[10 1000])
ylabel('ORN Gain (Hz/V)')
xlabel('Mean Light Power (\muW)')

% odour flicker light background -- LFP gain
subplot(3,3,2), hold on
x = []; y = []; base_gain = [];
for i = 1:length(paradigm)
	if light_background_paradigms(i)
		ci = find(paradigm(i) == plot_these_paradigms);
		plot(mean_light_power(i),LFP_gain(i),'k+');
		x = [x; mean_light_power(i)];
		y = [y; LFP_gain(i)];
		if ci == 1
			base_gain = [base_gain LFP_gain(i)];
		end
	end
end
plot([min(x) max(x)],[nanmean(base_gain) nanmean(base_gain)],'k--')
set(gca,'XScale','linear','YScale','log','YLim',[1e0 1e2])
ylabel('LFP Gain (mV/V)')
xlabel('Mean Light Power (\muW)')



% % ##       ####  ######   ##     ## ########    ########       ##  ######   
% % ##        ##  ##    ##  ##     ##    ##       ##            ##  ##    ##  
% % ##        ##  ##        ##     ##    ##       ##           ##   ##        
% % ##        ##  ##   #### #########    ##       ######      ##    ##   #### 
% % ##        ##  ##    ##  ##     ##    ##       ##         ##     ##    ##  
% % ##        ##  ##    ##  ##     ##    ##       ##        ##      ##    ##  
% % ######## ####  ######   ##     ##    ##       ##       ##        ######   

if ~exist('od','var')
	p = '/local-data/DA-paper/Chrimson/light-odor-combinations/light-flicker/v3';
	od = raw2ORNData(p,'led_power_func',light_power_fit,'use_led',true,'filter_LFP',false);

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
	od_odour = raw2ORNData(p,'led_power_func',light_power_fit,'use_led',false,'filter_LFP',false);
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

clear ax
ax(1) = subplot(3,3,5); hold on
ax(2) = subplot(3,3,6); hold on
c2 = lines(size(od,1));
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
			if i == 4
				axes(ax(1))
				[~,data] = plotPieceWiseLinear(x,y,'nbins',50,'Color',c(ci,:));
			end
			middle_segment = (y > max(y)/3 & y < 2*max(y)/3);
			ff = fit(x(middle_segment),y(middle_segment),'poly1');
			gain(j) = ff.p1;
			odour_levels(j) = nanmean((nanmean(od_odour(i,j).stimulus(:,:),2)));
			ci = ci+1;
		end
		rm_this = isnan(odour_levels) | isnan(gain) | gain == 0 ;
		odour_levels(rm_this) = []; gain(rm_this) = [];
		plot(ax(2),odour_levels,gain,'-+','Color',c2(i,:))
	end
end
xlabel(ax(1),'Projected Stimulus (\muW)')
ylabel(ax(1),'ab3A Firing Rate (Hz)')

set(ax(2),'XScale','log','YLim',[0 2],'XTick',[1e-4 1e-3 1e-2 1e-1 1e0])
xlabel(ax(2),'Odour Background (V)')
ylabel(ax(2),'ORN Gain (Hz/\muW)')

% add an inset showing how the odour affects the neuron
ax = axes(); hold on
ax.Position = [0.8 0.52 0.12 0.12];
ax.Tag = 'inset';
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
	plot(ax,s,f,'-+','Color',c(i,:))
end
set(ax,'XScale','log','XTick',[1e-4 1e-3 1e-2 1e-1 1e0])
ylabel(ax,'Firing rate (Hz)')

subplot(3,3,1), hold on
o = imread('../images/chrimson-1.png');
imagesc(o);
axis ij
axis image
axis off

subplot(3,3,4), hold on
o = imread('../images/chrimson-2.png');
imagesc(o);
axis ij
axis image
axis off

subplot(3,3,7), hold on
o = imread('../images/chrimson-3.png');
imagesc(o);
axis ij
axis image
axis off

prettyFig('fs',13,'FixLogX',true,'lw',1.3,'plw',1.2)
labelFigure


if being_published
	snapnow
	delete(gcf)
end


%% Version Info
%
pFooter;


%% Supplementary Figure
