% fig_light_mechanism.m
% 
% created by Srinivas Gorur-Shandilya at 9:50 , 28 December 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


%% Mechanism of Gain Control: Optogenetic Activation

pHeader;
dm = dataManager;

% build the LED map
x = [.75:0.05:1.1 1.2:.1:3.6 3.8 4 4.2 4.5];
x = [0:0.1:0.7 x];
y = [0 4 12.8 24.1 36.2 48.1 60 70 95 114 134 151 167 184 201 219 236 252 269 283 299 314 328 343 357 370 383 394 404 413 419 422 422 421 419 418 417];
y = [0*(0:0.1:0.7) y ];
light_power_fit = fit(x(:),y(:),'spline');

main_fig = figure('outerposition',[0 0 850 850],'PaperUnits','points','PaperSize',[850 850]); hold on
for i = 9:-1:1
	ax(i) = subplot(3,3,i); hold on
end

axs(5) = axes;
axs(5).Position = [0.8 0.53 0.1 0.1];
hold on

axs(3) = axes;
axs(3).Position = [0.8 0.85 0.1 0.1];
hold on

% we're killing the supplementary figure; replacing it with insets in the main figure
% supp_fig = figure('outerposition',[0 0 1100 700],'PaperUnits','points','PaperSize',[1100 700]); hold on
% for i = 6:-1:1
% 	axs(i) = subplot(2,3,i); hold on
% end

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

[~, ~, fA, paradigm, orn, fly, AllControlParadigms] = consolidateData(dm.getPath('38901017007c50ea52637891619ab91c'),true);


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

% also compute gains as ratio of std. devs. 
gain_s = NaN*gain;
for i = 1:length(gain)
	gain_s(i) = std(fA(a:z,i))/std(LED(a:z,i));
end
gain_s(gain_s == 0) = NaN;

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

% plot i/o curves for lowest and highest contrast
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
	%x = x - mean(x) + mean(mean(LED(a:z,light_s == contrast_levels(i))));

	[~,data] = plotPieceWiseLinear(x,y,'make_plot',false,'nbins',40);
	plot(ax(8),data.x,data.y,'Color',c(i,:),'LineWidth',3)
end
xlabel(ax(8),'Projected light stimulus (\muW)')
ylabel(ax(8),'ab3A firing rate (Hz)')
set(ax(8),'XLim',[0 500],'YLim',[0 50])
orn(find(m==2,1,'first')) = 6;

% plot gain vs. contrast
these_orns = unique(orn(m==2));
norns = length(these_orns); 
c = lines(norns);
for i = 1:length(c)
	plot(ax(9),light_s(m==2 & orn == these_orns(i)),gain(m==2 & orn == these_orns(i)),'-+','Color',c(i,:))
end
xlabel(ax(9),'\sigma_{Light} (\muW)')
ylabel(ax(9),'ab3A ORN gain (Hz/\muW)')
set(ax(9),'YLim',[0 .16],'XLim',[0 130])

% make the supp. figures. 
% for i = 1:length(c)
% 	plot(axs(6),light_s(m==2 & orn == these_orns(i)),gain(m==2 & orn == these_orns(i)),'-+','Color',c(i,:))
% end
% xlabel(axs(6),'\sigma_{Light} (\muW)')
% ylabel(axs(6),'\sigma_{Firing rate}/\sigma_{Stimulus} (Hz/\muW)')
% set(axs(6),'XLim',[30 130],'YLim',[0 .17])

%    ##       ####  ######   ##     ## ########    ########   ######   
%    ##        ##  ##    ##  ##     ##    ##       ##     ## ##    ##  
%    ##        ##  ##        ##     ##    ##       ##     ## ##        
%    ##        ##  ##   #### #########    ##       ########  ##   #### 
%    ##        ##  ##    ##  ##     ##    ##       ##     ## ##    ##  
%    ##        ##  ##    ##  ##     ##    ##       ##     ## ##    ##  
%    ######## ####  ######   ##     ##    ##       ########   ######   

clearvars -except ax axs being_published dm light_power_fit od od_odour uts main_fig supp_fig

[PID, LFP, fA, paradigm, orn, fly, AllControlParadigms, paradigm_hashes] = consolidateData(dm.getPath('9df40545ca5197f63f2dd4af7c905316'),true);
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


% extract filters and compute gain for LFP and firing rate
[K_LFP,LFP_pred,LFP_gain] = extractFilters(PID,filtered_LFP,'use_cache',true,'a',a,'z',z);
[K_fA,fA_pred,fA_gain] = extractFilters(PID,fA,'use_cache',true,'a',a,'z',z);

% also compute gains directly w/o LN model
LFP_gain_s = NaN*LFP_gain;
fA_gain_s = NaN*fA_gain;
for i = 1:length(fA_gain)
	LFP_gain_s(i) = std(LFP(a:z,i))/std(PID(a:z,i));
	fA_gain_s(i) = std(fA(a:z,i))/std(PID(a:z,i));
end
LFP_gain_s(LFP_gain_s == 0) = NaN;
fA_gain_s(fA_gain_s == 0) = NaN;


% plot the light background odour flicker cases -- total gain
plot_these_paradigms = sort(unique(paradigm(light_background_paradigms)));
these_orns = unique(orn(light_background_paradigms));
c = lines(length(these_orns));
for i = 1:length(these_orns)
	ii = orn == these_orns(i) & ismember(paradigm,plot_these_paradigms);
	if sum(ii) > 3
		plot(ax(3),mean_light_power(ii),fA_gain(ii),'-+','Color',c(i,:));
		%plot(axs(2),mean_light_power(ii),fA_gain_s(ii),'-+','Color',c(i,:));
	end
end

set(ax(3),'XScale','linear','YScale','log','YLim',[10 1000])
ylabel(ax(3),'ab3A ORN gain (Hz/V)')
xlabel(ax(3),'Mean light power (\muW)')
% set(axs(2),'XScale','linear','YScale','log','YLim',[10 1000])
% ylabel(axs(2),'\sigma_{Firing rate}/\sigma_{Stimulus} (Hz/V)')
% xlabel(axs(2),'Mean light power (\muW)')


% plot the transduction gain
for i = 1:length(these_orns)
	ii = orn == these_orns(i) & ismember(paradigm,plot_these_paradigms);
	if sum(ii) > 3
		plot(ax(2),mean_light_power(ii),LFP_gain(ii),'-+','Color',c(i,:));
		%plot(axs(1),mean_light_power(ii),10*LFP_gain_s(ii),'-+','Color',c(i,:));
	end
end

set(ax(2),'XScale','linear','YScale','log','YLim',[1 100])
ylabel(ax(2),'ab3 transduction gain (mV/V)')
xlabel(ax(2),'Mean light power (\muW)')
% set(axs(1),'XScale','linear','YScale','log','YLim',[1 100])
% ylabel(axs(1),'\sigma_{LFP}/\sigma_{Stimulus} (mV/V)')
% xlabel(axs(1),'Mean light power (\muW)')

% show that firing rate increases with increasing light stim
for i = 1:length(these_orns)
	ii = orn == these_orns(i) & ismember(paradigm,plot_these_paradigms);
	plot(axs(3),mean_light_power(ii),max(fA(1:5e3,ii)),'-+','Color',c(i,:));
end
ylabel(axs(3),['Light-induced' char(10) 'firing rate (Hz)'])
xlabel(axs(3),'Mean light power (\muW)')
set(axs(3),'XLim',[0 500],'YLim',[0 90])

% % ##       ####  ######   ##     ## ########    ########       ##  ######   
% % ##        ##  ##    ##  ##     ##    ##       ##            ##  ##    ##  
% % ##        ##  ##        ##     ##    ##       ##           ##   ##        
% % ##        ##  ##   #### #########    ##       ######      ##    ##   #### 
% % ##        ##  ##    ##  ##     ##    ##       ##         ##     ##    ##  
% % ##        ##  ##    ##  ##     ##    ##       ##        ##      ##    ##  
% % ######## ####  ######   ##     ##    ##       ##       ##        ######   

if ~exist('od','var')
	p = dm.getPath('6c02fb233f8221b2d22d98ebf8a48edb');
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
				axes(ax(5))
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
	end
	% weird outlier which is obviously wrong
	if i == 2 & length(gain) == 5
		gain(2) = []; odour_levels(2) = [];
	end
	plot(ax(6),odour_levels,gain,'-+','Color',c2(i,:))
end


xlabel(ax(5),'Projected Light Stimulus (\muW)')
ylabel(ax(5),'ab3A Firing Rate (Hz)')
set(ax(5),'XLim',[0 150],'YLim',[0 65])

set(ax(6),'XScale','log','YLim',[0 4],'XTick',[1e-4 1e-3 1e-2 1e-1 1e0],'XLim',[1e-4 1])
xlabel(ax(6),'Odor background (V)')
ylabel(ax(6),'ab3A ORN gain (Hz/\muW)')

% show that firing rate increases with odour in supp. figure. 
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
	plot(axs(5),s,f,'-+','Color',c(i,:))
end
set(axs(5),'XScale','log','XTick',[1e-4 1e-2 1e0])
ylabel(axs(5),['Odor-induced' char(10) 'firing rate (Hz)'])
xlabel(axs(5),'Odor Background (V)')

% % also show gain computed directly
% for i = 1:size(od,1)
% 	gain = NaN(size(od,2),1);
% 	odour_levels = NaN(size(od,2),1);
% 	for j = 1:size(od,2)
% 		if od(i,j).n_trials > 0
% 			x = (nanmean(od(i,j).firing_projected(uts,:),2)); x = x(:);
% 			y = (nanmean(od(i,j).firing_rate(uts,:),2)); y = y(:);
% 			gain(j) = std(y)/std(x);
% 			odour_levels(j) = nanmean((nanmean(od_odour(i,j).stimulus(:,:),2)));
% 		end
% 		rm_this = isnan(odour_levels) | isnan(gain) | gain == 0 ;
% 		odour_levels(rm_this) = []; gain(rm_this) = [];
% 		plot(axs(4),odour_levels,gain,'-+','Color',c2(i,:))
% 	end
% end
% set(axs(4),'XScale','log','XTick',[1e-4 1e-3 1e-2 1e-1 1e0],'YLim',[0 1],'XLim',[1e-4 1])
% ylabel(axs(4),'ab3A ORN gain (Hz/\muW)')
% xlabel(axs(4),'Odor background (V)')

axes(ax(1))
o = imread('../images/chrimson-1.png');
imagesc(o);
axis ij
axis image
axis off

axes(ax(4))
o = imread('../images/chrimson-2.png');
imagesc(o);
axis ij
axis image
axis off

axes(ax(7))
o = imread('../images/chrimson-3.png');
imagesc(o);
axis ij
axis image
axis off

% cosmetic fixes
ax(3).YLim = [50 5e3];

prettyFig(main_fig,'fs',13,'FixLogX',true,'lw',1.5,'plw',1.2)

labelFigure('ignore_these',[axs ax([1 4 7])])

return

deintersectAxes(ax([2 3 6]))
uistack(axs(3),'top')
uistack(axs(5),'top')

if being_published
	snapnow
	delete(main_fig)
end

% prettyFig(supp_fig,'fs',15,'FixLogX',true,'lw',1.7,'plw',2)

% if being_published
% 	snapnow
% 	delete(supp_fig)
% end


%% Version Info
%
pFooter;


%% Supplementary Figure
