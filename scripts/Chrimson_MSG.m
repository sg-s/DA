% Chrimson_MSG.m
% 
% created by Srinivas Gorur-Shandilya at 4:24 , 12 November 2015. Contact me at http://srinivas.gs/contact/
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

%% Gain changes with light in Chrimson flies
% In this document we look at how gain changes using changes in light mean and variance in ORNs expressing Chrimson. 


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


%% Gain changes with mean light. 
% In this section we look at how ORN gain changes with light background, where we measure the gain in response to a fluctuating light stimulus. No odour at all in this experiment. In the following figure, we plot input-output curves of the ORNs as we change the mean light intensity. 


% extract filters
[K,fp,gain,gain_err] = extractFilters(LED,fA,'use_cache',true,'a',a,'z',z);

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
nplots = length(unique(r));
all_r = unique(r);
for i = 1:nplots
	autoPlot(nplots,i); hold on
	plot_this = r == all_r(i);
	mean_stims = unique(m(plot_this));
	c = parula(length(mean_stims)+1);
	for j = 1:length(mean_stims)
		do_this = (plot_this & m == mean_stims(j));
		x = fp(:,do_this);
		x(:,nansum(x)==0) = [];
		y = fA(:,do_this);
		if ~isvector(x)
			x = mean2(x);
		end
		x = x + nanmean(nanmean(LED(a:z,do_this)));
		if ~isvector(y)
			y = mean2(y);
		end
		x = x(a:z); y = y(a:z);
		plotPieceWiseLinear(x,y,'nbins',30,'Color',c(j,:));
	end
	xlabel('Projected Light Power (\muW)')
	ylabel('ORN Response (Hz)')
	title(['\sigma = ' oval(mean(light_s(plot_this))) ' \muW^2'])
end
prettyFig()

if being_published
	snapnow
	delete(gcf)
end

%%
% We can clearly see that the I/O curves change. In the following plot, we compare how gain changes with the mean of the input: 

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
nplots = length(unique(r));
all_r = unique(r);
for i = 1:nplots
	autoPlot(nplots,i); hold on
	plot_this = r == all_r(i);
	mean_stim = light_m(plot_this);
	g = gain(plot_this);
	plot(mean_stim,g,'k+')
	title(['\sigma = ' oval(mean(light_s(plot_this))) ' \muW^2'])
	xlabel('Mean Light (\muW)')
	ylabel('Gain (Hz/\muW)')
	set(gca,'XScale','log','YScale','log')
end

prettyFig()

if being_published
	snapnow
	delete(gcf)
end

%% Gain Changes with Light contrast 
% Now we look at how gain changes with the contrast of the light. 

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(light_s(m==2),gain(m==2),'k+')
xlabel('\sigma_{light}(\muW)')
ylabel('Gain (Hz/\muW)')
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
[status,git_hash]=unix('git rev-parse HEAD');
if ~status
	disp(git_hash)
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
