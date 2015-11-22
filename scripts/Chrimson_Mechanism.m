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
	filtered_LFP(:,i) = bandPass(LFP(:,i),1000,10);
end


% extract filters everywhere for LFP and gain

[K_LFP,LFP_pred,LFP_gain] = extractFilters(PID,filtered_LFP,'use_cache',true,'a',a,'z',z);
[K_fA,fA_pred,fA_gain] = extractFilters(PID,fA,'use_cache',true,'a',a,'z',z);

% correct the fA_pred by the mean stimulus
for i = 1:length(paradigm)
	fA_pred(:,i) = fA_pred(:,i) + mean_PID(i);
end

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
x = []; y = [];
for i = 1:length(paradigm)
	if odour_background_paradigms(i)
		ci = find(paradigm(i) == plot_these_paradigms);
		plot(mean_PID(i),fA_gain(i),'+','Color',c(ci,:));
		x = [x mean_PID(i)];
		y = [y fA_gain(i)];
	end
end
set(gca,'XScale','log','YScale','log')
ff = fit(x(:),y(:),'power1','Lower',[-Inf -1],'Upper',[Inf -1]);
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
x = []; y = [];
for i = 1:length(paradigm)
	if odour_background_paradigms(i)
		ci = find(paradigm(i) == plot_these_paradigms);
		plot(mean_PID(i),LFP_gain(i),'+','Color',c(ci,:));
		x = [x mean_PID(i)];
		y = [y LFP_gain(i)];
	end
end
set(gca,'XScale','log','YScale','log')
ff = fit(x(:),y(:),'power1');
plot(sort(x),ff(sort(x)),'r')
ylabel('Gain (mV/V)')
xlabel('Mean Stimulus (V)')

prettyFig()

if being_published
	snapnow
	delete(gcf)
end

%% 
% So it looks like we can reproduce, quantitatively, the gain scaling properties we saw in WT ORNs in these ORNs that express Chrimson. 

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
x = []; y = [];
for i = 1:length(paradigm)
	if light_background_paradigms(i)
		ci = find(paradigm(i) == plot_these_paradigms);
		plot(mean_light_power(i),fA_gain(i),'+','Color',c(ci,:));
		x = [x mean_light_power(i)];
		y = [y fA_gain(i)];
	end
end
% set(gca,'XScale','log','YScale','log')
% x(x==0) = 1e-3;
% ff = fit(x(:),y(:),'power1','Lower',[-Inf -1],'Upper',[Inf -1]);
% plot(sort(x),ff(sort(x)),'r')
ylabel('Gain (Hz/V)')
xlabel('Mean Stimulus (V)')

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
xlabel('Projected Stimulus (V)')
ylabel('LFP Response (mV)')
subplot(2,2,4), hold on
x = []; y = [];
for i = 1:length(paradigm)
	if light_background_paradigms(i)
		ci = find(paradigm(i) == plot_these_paradigms);
		plot(mean_light_power(i),LFP_gain(i),'+','Color',c(ci,:));
		x = [x mean_light_power(i)];
		y = [y LFP_gain(i)];
	end
end
% set(gca,'XScale','log','YScale','log')
% x(x==0) = 1e-3;
% ff = fit(x(:),y(:),'power1');
% plot(sort(x),ff(sort(x)),'r')
ylabel('Gain (mV/V)')
xlabel('Mean Stimulus (\muW)')

prettyFig()

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

t = toc;

%% 
% This document was built in: 
disp(strcat(oval(t,3),' seconds.'))

% tag the file as being published 

if being_published
	unix(['tag -a published ',which(mfilename)]);
	unix(['tag -r publish-failed ',which(mfilename)]);
end
