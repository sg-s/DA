% fig_light_mechanism.m
% 
% created by Srinivas Gorur-Shandilya at 9:50 , 28 December 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


%% Mechanism of Gain Control: Optogenetic Activation




% add homebrew path
path1 = getenv('PATH');
if isempty(strfind(path1,':/usr/local/bin'))
    path1 = [path1 ':/usr/local/bin'];
end
setenv('PATH', path1);

% this code determines if this function is being called
calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
		unix(['tag -a publish-failed ',which(mfilename)]);
		unix(['tag -r published ',which(mfilename)]);
	end
end
tic


figure('outerposition',[0 0 1000 900],'PaperUnits','points','PaperSize',[1000 900]); hold on

%     ##       ####  ######   ##     ## ########    ##     ##  ######   ######   
%     ##        ##  ##    ##  ##     ##    ##       ###   ### ##    ## ##    ##  
%     ##        ##  ##        ##     ##    ##       #### #### ##       ##        
%     ##        ##  ##   #### #########    ##       ## ### ##  ######  ##   #### 
%     ##        ##  ##    ##  ##     ##    ##       ##     ##       ## ##    ##  
%     ##        ##  ##    ##  ##     ##    ##       ##     ## ##    ## ##    ##  
%     ######## ####  ######   ##     ##    ##       ##     ##  ######   ######   


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


subplot(3,3,5), hold on
plot_this = r == 1;
mean_stim = light_m(plot_this);
g = gain(plot_this);
plot(mean_stim,g,'k+')
xlabel('Mean Light Power(\muW)')
ylabel('ORN Gain (Hz/\muW)')
set(gca,'XScale','log','YScale','log')
mean_stim = mean_stim(~isnan(g));
g = g(~isnan(g));
ff = fit(mean_stim(:),g(:),'power1');%,'Lower',[-Inf -1],'Upper',[Inf -1]);
plot([min(mean_stim) max(mean_stim)],ff([min(mean_stim) max(mean_stim)]),'k--')
ff = fit(mean_stim(:),g(:),'power1','Lower',[-Inf -1],'Upper',[Inf -1]);
plot([min(mean_stim) max(mean_stim)],ff([min(mean_stim) max(mean_stim)]),'r')

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


subplot(3,3,7); hold on
ylabel('Probability')
xlabel('Light Power (\muW)')

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
	plot(xx(2:end),yy,'Color',c(i,:))
end

subplot(3,3,8), hold on
plot(light_s(m==2),gain(m==2),'k+')
xlabel('\sigma_{light}(\muW)')
ylabel('ORN Gain (Hz/\muW)')

%  #######  ########   #######  ##     ## ########  
% ##     ## ##     ## ##     ## ##     ## ##     ## 
% ##     ## ##     ## ##     ## ##     ## ##     ## 
% ##     ## ##     ## ##     ## ##     ## ########  
% ##     ## ##     ## ##     ## ##     ## ##   ##   
% ##     ## ##     ## ##     ## ##     ## ##    ##  
%  #######  ########   #######   #######  ##     ## 

%  #######  ########   #######  ##     ## ########  
% ##     ## ##     ## ##     ## ##     ## ##     ## 
% ##     ## ##     ## ##     ## ##     ## ##     ## 
% ##     ## ##     ## ##     ## ##     ## ########  
% ##     ## ##     ## ##     ## ##     ## ##   ##   
% ##     ## ##     ## ##     ## ##     ## ##    ##  
%  #######  ########   #######   #######  ##     ## 

%  ######   #######  ##    ## ######## ########   #######  ##        ######  
% ##    ## ##     ## ###   ##    ##    ##     ## ##     ## ##       ##    ## 
% ##       ##     ## ####  ##    ##    ##     ## ##     ## ##       ##       
% ##       ##     ## ## ## ##    ##    ########  ##     ## ##        ######  
% ##       ##     ## ##  ####    ##    ##   ##   ##     ## ##             ## 
% ##    ## ##     ## ##   ###    ##    ##    ##  ##     ## ##       ##    ## 
%  ######   #######  ##    ##    ##    ##     ##  #######  ########  ######  



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
		filtered_LFP(:,i) = bandPass(LFP(:,i),1000,10);
	catch
	end
end


% extract filters everywhere for LFP and gain

[K_LFP,LFP_pred,LFP_gain] = extractFilters(PID,filtered_LFP,'use_cache',true,'a',a,'z',z);
[K_fA,fA_pred,fA_gain] = extractFilters(PID,fA,'use_cache',true,'a',a,'z',z);


% plot the firing gain for odour and odour
subplot(3,3,1), hold on
plot_these_paradigms = sort(unique(paradigm(odour_background_paradigms)));
x = []; y = []; base_gain = [];
for i = 1:length(paradigm)
	if odour_background_paradigms(i)
		ci = find(paradigm(i) == plot_these_paradigms);
		plot(mean_PID(i),fA_gain(i),'k+');
		x = [x mean_PID(i)];
		y = [y fA_gain(i)];
		if ci == 1
			base_gain = [base_gain fA_gain(i)];
		end
	end
end
set(gca,'XScale','log','YScale','log')
rm_this = isnan(x) | isnan(y);
x(rm_this) = []; y(rm_this) = [];
ff = fit(x(:),y(:),'power1','Lower',[-Inf -1],'Upper',[Inf -1]);
% also plot a flat line for comparison
plot([min(x) max(x)],[nanmean(base_gain) nanmean(base_gain)],'k--')
plot(sort(x),ff(sort(x)),'r')
ylabel('ORN Gain (Hz/V)')
xlabel('Odor concentration (V)')


% also plot the LFP gain
subplot(3,3,3), hold on
x = []; y = []; base_gain = [];
for i = 1:length(paradigm)
	if odour_background_paradigms(i)
		ci = find(paradigm(i) == plot_these_paradigms);
		plot(mean_PID(i),LFP_gain(i),'k+');
		x = [x mean_PID(i)];
		y = [y LFP_gain(i)];
		if ci == 1
			base_gain = [base_gain LFP_gain(i)];
		end
	end
end
plot([min(x) max(x)],[nanmean(base_gain) nanmean(base_gain)],'k--')
set(gca,'XScale','log','YScale','log','YLim',[.05 2])
rm_this = isnan(x) | isnan(y);
x(rm_this) = []; y(rm_this) = [];
ff = fit(x(:),y(:),'power1','Lower',[-Inf -1],'Upper',[Inf -1]);
plot(sort(x),ff(sort(x)),'r')
ylabel('LFP Gain (mV/V)')
xlabel('Odor concentration (V)')

%    ##       ####  ######   ##     ## ########    ########   ######   
%    ##        ##  ##    ##  ##     ##    ##       ##     ## ##    ##  
%    ##        ##  ##        ##     ##    ##       ##     ## ##        
%    ##        ##  ##   #### #########    ##       ########  ##   #### 
%    ##        ##  ##    ##  ##     ##    ##       ##     ## ##    ##  
%    ##        ##  ##    ##  ##     ##    ##       ##     ## ##    ##  
%    ######## ####  ######   ##     ##    ##       ########   ######   

subplot(3,3,4), hold on
plot_these_paradigms = sort(unique(paradigm(light_background_paradigms)));

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
% x(x==0) = 1e-3; x = x(:);
% ff = fit(x(~isnan(y)),y(~isnan(y)),'power1','Lower',[-Inf -1],'Upper',[Inf -1]);
% plot(sort(x),ff(sort(x)),'r')
ylabel('ORN Gain (Hz/V)')
xlabel('Mean Light Power (\muW)')


subplot(3,3,6), hold on
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
set(gca,'XScale','linear','YScale','log','YLim',[1e-2 1e1])
% x(x==0) = 1e-3;
% ff = fit(x(:),y(:),'power1');
% plot(sort(x),ff(sort(x)),'r')
ylabel('LFP Gain (mV/V)')
xlabel('Mean Light Power (\muW)')


% ##       ####  ######   ##     ## ########    ########       ##  ######   
% ##        ##  ##    ##  ##     ##    ##       ##            ##  ##    ##  
% ##        ##  ##        ##     ##    ##       ##           ##   ##        
% ##        ##  ##   #### #########    ##       ######      ##    ##   #### 
% ##        ##  ##    ##  ##     ##    ##       ##         ##     ##    ##  
% ##        ##  ##    ##  ##     ##    ##       ##        ##      ##    ##  
% ######## ####  ######   ##     ##    ##       ##       ##        ######   


p = '/local-data/DA-paper/Chrimson/light-odor-combinations/light-flicker/ok/';
[PID, ~, fA, paradigm, orn, fly, AllControlParadigms] = consolidateData(p,true);
a = 15e3; z = 45e3;

rm_this = isnan(sum(fA(a:z,:))) | ~nansum(fA(a:z,:));
fA(:,rm_this) = [];
PID(:,rm_this) = [];
paradigm(rm_this) = [];
orn(rm_this) = [];
fly(rm_this) = [];
LED = 0*PID;

% make a vector showing which odour we use, and what the concentration is. 
% also make a LED stimulus vector
odour = {}; conc = 0*paradigm; nominal_conc = 0*paradigm;
for i = 1:length(paradigm)
	if any(strfind(AllControlParadigms(paradigm(i)).Name,'back'))
		odour{i} = '2ac';
	end
	if any(strfind(AllControlParadigms(paradigm(i)).Name,'2ac'))
		odour{i} = '2ac';
	end
	if any(strfind(AllControlParadigms(paradigm(i)).Name,'ethanol'))
		odour{i} = 'ethanol';
	end
	if str2double(AllControlParadigms(paradigm(i)).Name(end)) > 0
		conc(i) = nanmean(PID(a:z,i));
		nominal_conc(i)  =str2double(AllControlParadigms(paradigm(i)).Name(strfind(AllControlParadigms(paradigm(i)).Name,'=')+1:end));
	end
	% figure out which output is the LED  -- probably the one with the highest variance
	[~,load_this] = max(std(AllControlParadigms(paradigm(i)).Outputs'));
	LED(:,i) = light_power_fit(AllControlParadigms(paradigm(i)).Outputs(load_this,1:10:end));
end

[K_fA,fA_pred,fA_gain] = extractFilters(LED,fA,'use_cache',true,'a',a,'z',z);
odour_levels = sort(unique(nominal_conc(find(strcmp(odour,'2ac')))));
c = parula(length(odour_levels)+1);

subplot(3,3,2), hold on

y = fA_gain(nominal_conc == 0);
x = 0*y;
plot(x,y,'k+';
for i = 2:length(odour_levels)
	x = mean(PID(a:z,intersect(find(nominal_conc == odour_levels(i)),find(strcmp(odour,'2ac')))));
	y = (fA_gain(intersect(find(nominal_conc == odour_levels(i)),find(strcmp(odour,'2ac')))));
	plot(x,y,'k+');
end
xlabel('Odor concentration (V)')
ylabel('ORN Gain (Hz/\muW)')
set(gca,'YLim',[.05 .2])


prettyFig('fs=15;','FixLogX=true;')

if being_published
	snapnow
	delete(gcf)
end


%% Version Info
% The file that generated this document is called:
disp(mfilename)

%%
% and its md5 hash is:
Opt.Input = 'file';
disp(dataHash(strcat(mfilename,'.m'),Opt))

%%
% This file should be in this commit:
[status,m]=unix('git rev-parse HEAD');
if ~status
	disp(m)
end

%%
% This file has the following external dependencies:
showDependencyHash(mfilename);

t = toc;

%% 
% This document was built in: 
disp(strcat(oval(t,3),' seconds.'))

% tag the file as being published 

if being_published
	unix(['tag -a published ',which(mfilename)]);
	unix(['tag -r publish-failed ',which(mfilename)]);
end
