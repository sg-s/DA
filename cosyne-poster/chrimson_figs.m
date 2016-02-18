% chrimson_figs.m
% this script makes figures for the cosyne poster, showing gain control by light activation
% 
% created by Srinivas Gorur-Shandilya at 9:50 , 28 December 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


%% Mechanism of Gain Control: Optogenetic Activation


% make figure placeholders 
fig_handle=figure('outerposition',[0 0 1200 600],'PaperUnits','points','PaperSize',[1200 600]); hold on

for i = 1:8
	axes_handles(i) = subplot(2,4,i);
	hold(axes_handles(i),'on');
end
axes_handles(5) = subplot(2,2,3);
hold(axes_handles(5),'on');


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


x = [.75:0.05:1.1 1.2:.1:3.6 3.8 4 4.2 4.5];
x = [0:0.1:0.7 x];
y = [0 4 12.8 24.1 36.2 48.1 60 70 95 114 134 151 167 184 201 219 236 252 269 283 299 314 328 343 357 370 383 394 404 413 419 422 422 421 419 418 417];
y = [0*(0:0.1:0.7) y ];
light_power_fit = fit(x(:),y(:),'spline');

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
time = 1e-3*(1:length(LED));

% plot the light stim
contrast_levels = sort(unique(light_s(m==2)));
% skip the third
contrast_levels(3) = [];
c = parula(3);
c(3,:) = [218 142 0]/255; % yellow too bright

for i = 1:length(contrast_levels)
	plot_these = light_s == contrast_levels(i);
	y = nanmean(LED(a:z,plot_these),2);
	y = y(1:5e3);
	y(1:1e3) = NaN;
	x = 25 + 1e-3*(1:5e3) + i*5;
	plot(axes_handles(5),x,y,'Color',c(i,:))
end

set(axes_handles(5),'XTick',[31:2:35 36:2:40 41:2:45 46:2:50],'XTickLabel',{'26','28','30','26','28','30','26','28','30','26','28','30'})

% plot i/o curves

light_m = NaN*contrast_levels;

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

	[~,data] = plotPieceWiseLinear(x,y,'make_plot',false,'nbins',40);
	plot(axes_handles(7),data.x,data.y,'Color',c(i,:),'LineWidth',3)
end

% plot gain vs. contrast
plot(axes_handles(8),light_s(m==2),gain(m==2),'k+')



%    ##       ####  ######   ##     ## ########    ########   ######   
%    ##        ##  ##    ##  ##     ##    ##       ##     ## ##    ##  
%    ##        ##  ##        ##     ##    ##       ##     ## ##        
%    ##        ##  ##   #### #########    ##       ########  ##   #### 
%    ##        ##  ##    ##  ##     ##    ##       ##     ## ##    ##  
%    ##        ##  ##    ##  ##     ##    ##       ##     ## ##    ##  
%    ######## ####  ######   ##     ##    ##       ########   ######   

clearvars -except axes_handles time light_power_fit

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
for i = 1:width(LFP)
	try
		% band pass all the LFP
		filtered_LFP(:,i) = filtered_LFP(:,i) - fastFiltFilt(ones(1e4,1),1e4,filtered_LFP(:,i));
		filtered_LFP(:,i) = filtered_LFP(:,i)*10; % to get the units right, now in mV
	catch
	end
end


% extract filters everywhere for LFP and gain
[K_LFP,LFP_pred,LFP_gain] = extractFilters(PID,filtered_LFP,'use_cache',true,'a',a,'z',z);
[K_fA,fA_pred,fA_gain] = extractFilters(PID,fA,'use_cache',true,'a',a,'z',z);

% plot the odour stimulus
plot(axes_handles(1),time(4e3:12e3),nanmean(PID(4e3:12e3,light_background_paradigms),2),'LineWidth',3,'Color','k')

% plot the light stimulus
light_bg_levels = unique(mean_light_power(light_background_paradigms));
c = parula(length(light_bg_levels)+1);
for i = 1:3:length(light_bg_levels)
	y = light_bg_levels(i);
	y = y*ones(8e3+1,1);
	y(1:1e3) = 0;
	plot(axes_handles(2),time(4e3:12e3),y,'Color',c(i,:))
end


% plot the firing gain for odour with light bg
x = []; y = []; base_gain = [];
for i = 1:length(paradigm)
	if light_background_paradigms(i)
		plot(axes_handles(4),mean_light_power(i),fA_gain(i),'+','Color',c(find(light_bg_levels == mean_light_power(i)),:),'LineWidth',2);
		x = [x mean_light_power(i)];
		y = [y fA_gain(i)];
		if mean_light_power(i) == 0
			base_gain = [base_gain fA_gain(i)];
		end
	end
end
set(axes_handles(4),'XScale','linear','YScale','log')
% also plot a flat line for comparison
plot(axes_handles(4),[min(x) max(x)],[nanmean(base_gain) nanmean(base_gain)],'k--')


% plot the LFP gain for odour with light bg
x = []; y = []; base_gain = [];
for i = 1:length(paradigm)
	if light_background_paradigms(i)
		plot(axes_handles(3),mean_light_power(i),LFP_gain(i),'+','Color',c(find(light_bg_levels == mean_light_power(i)),:),'LineWidth',2);
		x = [x mean_light_power(i)];
		y = [y LFP_gain(i)];
		if mean_light_power(i) == 0
			base_gain = [base_gain LFP_gain(i)];
		end
	end
end
set(axes_handles(3),'XScale','linear','YScale','log')
% also plot a flat line for comparison
plot(axes_handles(3),[min(x) max(x)],[nanmean(base_gain) nanmean(base_gain)],'k--')

% CREATE placeholders for all x and y labels
for i = 1:length(axes_handles)
	try
		xlabel(axes_handles(i),'A A')
		ylabel(axes_handles(i),'A A')
	end
end

prettyFig('fs=18;','FixLogX=true;','lw=2;','plw=2;')
