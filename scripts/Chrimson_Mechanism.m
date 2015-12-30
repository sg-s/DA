% Chrimson_Mechanism.m
% 
% created by Srinivas Gorur-Shandilya at 6:43 , 20 November 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


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

%% Light Stimulus
% First we see how the applied driver voltage to the LED is transformed into actual light power as measured by a light meter. 

x = [.75:0.05:1.1 1.2:.1:3.6 3.8 4 4.2 4.5];
x = [0:0.1:0.7 x];
y = [0 4 12.8 24.1 36.2 48.1 60 70 95 114 134 151 167 184 201 219 236 252 269 283 299 314 328 343 357 370 383 394 404 413 419 422 422 421 419 418 417];
y = [0*(0:0.1:0.7) y ];
light_power_fit = fit(x(:),y(:),'spline');

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
xx = 0:1e-2:5;
plot(xx,light_power_fit(xx),'r')
plot(x,y,'k+')
xlabel('Driving Voltage (V)')
ylabel('Power @ 627nm (\muW)')
prettyFig()

if being_published
	snapnow
	delete(gcf)
end


%% Odour Flicker with Odour Backgrounds
% In this section we study how the gain changes with background odour, to a odour flicker, in 22a-GAL4/+;UAS-Chrimson/+ flies. This serves as a control, and checks if these flies behave normally with the odour. 

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

% correct the fA_pred by the mean stimulus
for i = 1:length(paradigm)
	fA_pred(:,i) = fA_pred(:,i) + mean_PID(i);
end


%  #######  ########   #######  ########           #######  ########   #######  ########  
% ##     ## ##     ## ##     ## ##     ##         ##     ## ##     ## ##     ## ##     ## 
% ##     ## ##     ## ##     ## ##     ##         ##     ## ##     ## ##     ## ##     ## 
% ##     ## ##     ## ##     ## ########  ####### ##     ## ##     ## ##     ## ########  
% ##     ## ##     ## ##     ## ##   ##           ##     ## ##     ## ##     ## ##   ##   
% ##     ## ##     ## ##     ## ##    ##          ##     ## ##     ## ##     ## ##    ##  
%  #######  ########   #######  ##     ##          #######  ########   #######  ##     ## 


figure('outerposition',[0 0 700 700],'PaperUnits','points','PaperSize',[700 700]); hold on
subplot(2,2,1), hold on
plot_these_paradigms = sort(unique(paradigm(odour_background_paradigms)));
c = parula(length(plot_these_paradigms)+1);
for i = 1:length(plot_these_paradigms)
	x = fA_pred(a:z,paradigm == plot_these_paradigms(i));
	y = fA(a:z,paradigm == plot_these_paradigms(i));
	if ~isvector(x)
		x = mean2(x);
	end
	if ~isvector(y)
		y = mean2(y);
	end
	plotPieceWiseLinear(x,y,'nbins',30,'Color',c(i,:));
end
xlabel('Projected Stimulus (V)')
ylabel('ORN Response (Hz)')
subplot(2,2,2), hold on
x = []; y = []; base_gain = [];
for i = 1:length(paradigm)
	if odour_background_paradigms(i)
		ci = find(paradigm(i) == plot_these_paradigms);
		plot(mean_PID(i),fA_gain(i),'+','Color',c(ci,:));
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
ylabel('Gain (Hz/V)')
xlabel('Mean Stimulus (V)')

subplot(2,2,3), hold on
plot_these_paradigms = sort(unique(paradigm(odour_background_paradigms)));
c = parula(length(plot_these_paradigms)+1);
for i = 1:length(plot_these_paradigms)
	x = LFP_pred(a:z,paradigm == plot_these_paradigms(i));
	y = filtered_LFP(a:z,paradigm == plot_these_paradigms(i));
	if ~isvector(x)
		x = mean2(x);
	end
	if ~isvector(y)
		y = mean2(y);
	end
	plotPieceWiseLinear(x,y,'nbins',30,'Color',c(i,:));
end
xlabel('Projected Stimulus (V)')
ylabel('LFP Response (mV)')
subplot(2,2,4), hold on
x = []; y = []; base_gain = [];
for i = 1:length(paradigm)
	if odour_background_paradigms(i)
		ci = find(paradigm(i) == plot_these_paradigms);
		plot(mean_PID(i),LFP_gain(i),'+','Color',c(ci,:));
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
ff = fit(x(:),y(:),'power1');
plot(sort(x),ff(sort(x)),'r')
ylabel('Gain (mV/V)')
xlabel('Mean Stimulus (V)')
suptitle('Ethyl Acetate foreground and background in w;22a-GAL4/+;UAS-Chrimson/+ flies')
prettyFig('fs=12;')

if being_published
	snapnow
	delete(gcf)
end


%% 
% So it looks like we can reproduce, quantitatively, the gain scaling properties we saw in WT ORNs in these ORNs that express Chrimson. 

% ##       ####  ######   ##     ## ######## 
% ##        ##  ##    ##  ##     ##    ##    
% ##        ##  ##        ##     ##    ##    
% ##        ##  ##   #### #########    ##    
% ##        ##  ##    ##  ##     ##    ##    
% ##        ##  ##    ##  ##     ##    ##    
% ######## ####  ######   ##     ##    ##    

% ########   ######   
% ##     ## ##    ##  
% ##     ## ##        
% ########  ##   #### 
% ##     ## ##    ##  
% ##     ## ##    ##  
% ########   ######   



%% Gain Scaling with background light
% Now that we have established that these ORNs can change gain with supplemental odour, just like the WT ORNs, we see how they change gain to a fluctuating odour stimulus with supplemental light drive. 

figure('outerposition',[0 0 700 700],'PaperUnits','points','PaperSize',[700 700]); hold on
subplot(2,2,1), hold on
plot_these_paradigms = sort(unique(paradigm(light_background_paradigms)));
c = parula(length(plot_these_paradigms)+1);
for i = 1:length(plot_these_paradigms)
	x = fA_pred(a:z,paradigm == plot_these_paradigms(i));
	y = fA(a:z,paradigm == plot_these_paradigms(i));
	if ~isvector(x)
		x = mean2(x);
	end
	if ~isvector(y)
		y = mean2(y);
	end
	plotPieceWiseLinear(x,y,'nbins',30,'Color',c(i,:));
end
xlabel('Projected Stimulus (V)')
ylabel('ORN Response (Hz)')
subplot(2,2,2), hold on
x = []; y = []; base_gain = [];
for i = 1:length(paradigm)
	if light_background_paradigms(i)
		ci = find(paradigm(i) == plot_these_paradigms);
		plot(mean_light_power(i),fA_gain(i),'+','Color',c(ci,:));
		x = [x mean_light_power(i)];
		y = [y fA_gain(i)];
		if ci == 1
			base_gain = [base_gain fA_gain(i)];
		end
	end
end
plot([min(x) max(x)],[nanmean(base_gain) nanmean(base_gain)],'k--')
set(gca,'XScale','linear','YScale','log','YLim',[10 1000])
% x(x==0) = 1e-3;
% ff = fit(x(:),y(:),'power1','Lower',[-Inf -1],'Upper',[Inf -1]);
% plot(sort(x),ff(sort(x)),'r')
ylabel('ORN Gain (Hz/V)')
xlabel('Mean Light Power (\muW)')

subplot(2,2,3), hold on
plot_these_paradigms = sort(unique(paradigm(light_background_paradigms)));
c = parula(length(plot_these_paradigms)+1);
for i = 1:length(plot_these_paradigms)
	x = LFP_pred(a:z,paradigm == plot_these_paradigms(i));
	y = filtered_LFP(a:z,paradigm == plot_these_paradigms(i));
	if ~isvector(x)
		x = mean2(x);
	end
	if ~isvector(y)
		y = mean2(y);
	end
	plotPieceWiseLinear(x,y,'nbins',30,'Color',c(i,:));
end
xlabel('Projected Odor Stimulus (V)')
xlabel('Mean Light Power (\muW)')
subplot(2,2,4), hold on
x = []; y = []; base_gain = [];
for i = 1:length(paradigm)
	if light_background_paradigms(i)
		ci = find(paradigm(i) == plot_these_paradigms);
		plot(mean_light_power(i),LFP_gain(i),'+','Color',c(ci,:));
		x = [x mean_light_power(i)];
		y = [y LFP_gain(i)];
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
suptitle('Ethyl Acetate foreground, light background in w;22a-GAL4/+;UAS-Chrimson/+ flies')
prettyFig('fs=12;')

if being_published
	snapnow
	delete(gcf)
end

%%
% We don't see a similar gain scaling. Perhaps the background light isn't strong enough? But this range of background light was sufficient to induce gain control in the same ORNs when the foreground fluctuation was also light. We are also driving the ORNs so hard that they pinch, so the background light is probably quite strong. 

%% 
% One possibility is that these ORNs adapt to the background light completely, so in effect they don't see it, but this is hard to reconcile with the fact that the same ORNs show Weber-like gain control with light foreground and light background. 

%% Light Flicker with Odour backgrounds.
% We now do the complementary experiment where we have a light flicker with an odour background, and see if the gain changes with odour background. In the first round of experiments, we used ethyl acetate as the background odour. The following figure shows how the input-output curve of the ORN (calcualted vs. light) varies with the background ethyl acetate odourant. 

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

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
x = mean(mean(LED(a:z,:))) + mean2(fA_pred(a:z,nominal_conc == 0));
y = mean2(fA(a:z,nominal_conc == 0));
plotPieceWiseLinear(x,y,'nbins',30,'Color',c(1,:));
xlabel('Projected Stimulus (\muW)')
ylabel('ORN Response (Hz)')
for i = 2:length(odour_levels)
	x = (fA_pred(a:z,intersect(find(nominal_conc == odour_levels(i)),find(strcmp(odour,'2ac')))));
	y = (fA(a:z,intersect(find(nominal_conc == odour_levels(i)),find(strcmp(odour,'2ac')))));
	if ~isvector(x)
		x = mean2(x);
	end
	if ~isvector(y)
		y = mean2(y);
	end
	plotPieceWiseLinear( mean(mean(LED(a:z,:))) + x,y,'nbins',30,'Color',c(i,:));
end

subplot(1,2,2), hold on
y = fA_gain(nominal_conc == 0);
x = 0*y;
plot(x,y,'+','Color',c(1,:));
for i = 2:length(odour_levels)
	x = mean(PID(a:z,intersect(find(nominal_conc == odour_levels(i)),find(strcmp(odour,'2ac')))));
	y = (fA_gain(intersect(find(nominal_conc == odour_levels(i)),find(strcmp(odour,'2ac')))));
	plot(x,y,'+','Color',c(i,:));
end
xlabel('Ethyl acetate concentration (V)')
ylabel('Gain (Hz/\muW)')
set(gca,'YLim',[.05 .2])

prettyFig('fs=14;')

if being_published
	snapnow
	delete(gcf)
end

%%
% A problem with using ethyl acetate was that the ORN was very sensitive to it, and pinched rapidly. To work around this, we instead used ethanol, which also excites the same neuron, but less effectively. This allowed us a finer control on the background activity of the neuron, and seemed to prevent pinching. 

odour_levels = sort(unique(nominal_conc(find(strcmp(odour,'ethanol')))));
c = parula(length(odour_levels)+1);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
x = mean(mean(LED(a:z,:))) + mean2(fA_pred(a:z,nominal_conc == 0));
y = mean2(fA(a:z,nominal_conc == 0));
plotPieceWiseLinear(x,y,'nbins',30,'Color',c(1,:));
xlabel('Projected Stimulus (\muW)')
ylabel('ORN Response (Hz)')
for i = 2:length(odour_levels)
	x = (fA_pred(a:z,intersect(find(nominal_conc == odour_levels(i)),find(strcmp(odour,'ethanol')))));
	y = (fA(a:z,intersect(find(nominal_conc == odour_levels(i)),find(strcmp(odour,'ethanol')))));
	if ~isvector(x)
		x = mean2(x);
	end
	if ~isvector(y)
		y = mean2(y);
	end
	plotPieceWiseLinear( mean(mean(LED(a:z,:))) + x,y,'nbins',30,'Color',c(i,:));
end

subplot(1,2,2), hold on
y = fA_gain(nominal_conc == 0);
x = 0*y;
plot(x,y,'+','Color',c(1,:));
for i = 2:length(odour_levels)
	x = mean(PID(a:z,intersect(find(nominal_conc == odour_levels(i)),find(strcmp(odour,'ethanol')))));
	y = (fA_gain(intersect(find(nominal_conc == odour_levels(i)),find(strcmp(odour,'ethanol')))));
	plot(x,y,'+','Color',c(i,:));
end
xlabel('Ethanol concentration (V)')
ylabel('Gain (Hz/\muW)')
set(gca,'YLim',[.05 .2])

prettyFig('fs=14;')

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
