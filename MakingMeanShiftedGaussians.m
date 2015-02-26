% MakingMeanShiftedGaussians.m
% in this document, we investigate methods to make mean shifted gaussians. 
% 
% created by Srinivas Gorur-Shandilya at 9:39 , 26 February 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end


%% Section 1: Presenting uniformly distributed noise 
% In this section, we drive the MFC with some uniformly distributed noise with a 100ms correlation time, being careful to choose the noise in the dilution space. 

load('/local-data/DA-paper/mean-shifted-gaussians/2015_02_26_odour_test_1.mat')
time = 1e-4*(1:length(data(2).PID));

figure('outerposition',[0 0 1400 1000],'PaperUnits','points','PaperSize',[1400 1000]); hold on


subplot(2,9,1:5), hold on
plot(time,data(2).MFC500')
set(gca,'XLim',[0 60])
xlabel('Time (s)')
ylabel('MFC Flow V)')

subplot(2,9,6:7), hold on
[r2,s] = rsquare(data(2).MFC500);
imagescnan(r2)
caxis([0 1])
colorbar
axis image
axis off
title(strcat('mean r^2 = ',oval(mean(r2(~isnan(r2))),2)))

subplot(2,9,8:9), hold on
imagescnan(s)
colorbar
axis image
axis off
title(strcat('mean slope = ',oval(mean(s(~isnan(s))),2)))


subplot(2,9,10:14), hold on
plot(time,data(2).PID')
set(gca,'XLim',[0 60])
xlabel('Time (s)')
ylabel('PID (V)')

subplot(2,9,15:16), hold on
[r2,s] = rsquare(data(2).PID);
imagescnan(r2)
caxis([0 1])
axis image
colorbar
axis off
title(strcat('mean r^2 = ',oval(mean(r2(~isnan(r2))),2)))

subplot(2,9,17:18), hold on
imagescnan(s)
colorbar
axis image
axis off
title(strcat('mean slope = ',oval(mean(s(~isnan(s))),2)))


PrettyFig;

if being_published
	snapnow
	delete(gcf)
end


%%
% First off, there is a long initial transient in the PID for the first trial, but subsequent trials are fine. Armed with this knowledge, we ignore the first trial. 





%% Version Info
% The file that generated this document is called:
disp(mfilename)

%%
% and its md5 hash is:
Opt.Input = 'file';
disp(DataHash(strcat(mfilename,'.m'),Opt))

%%
% This file should be in this commit:
[status,m]=unix('git rev-parse HEAD');
if ~status
	disp(m)
end