% MakingMSG.m
% documents attempts to construct gaussian distributed stimuli
% 
% created by Srinivas Gorur-Shandilya at 2:30 , 21 July 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


% this code determines if this function is being called by publish() or not
calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end
tic

%% Naive Assumptions
% First, we assume that the PID simply measures the concentration of the odour in gas phase, and that this simply depends on the fractional flow through the odourised bottle. In the following figures, we use ethyl acetate in a 30mL scintillation vial, with a 500mL/min MFC driving air through it. Main air is constant at 2L/min. The total flow is variable, as the odour flow fluctuates. 

%%
% In the following figure, we plot the distributions of MFC flow rates and the corresponding stimulus measurements:

load('/local-data/DA-paper/MSG-construction/2015_07_21_MSG_test_new_odour.mat')
haz_data = findData(data);

% trim the data
for i = haz_data
	data(i).PID = data(i).PID(:,20e4:10:55e4)';
	data(i).MFC500 = data(i).MFC500(:,20e4:10:55e4)';
end

c = parula(length(haz_data)+1);
figure('outerposition',[0 0 700 700],'PaperUnits','points','PaperSize',[700 700]); hold on
subplot(2,1,1), hold on
cc = 1;
for i = haz_data
	hist_this = data(i).MFC500;
	xx =  linspace(min(min(hist_this)),max(max(hist_this)),50);
	y = NaN(width(hist_this),50);
	for j = 1:width(hist_this)
		y(j,:) = hist(hist_this(:,j),xx);
		y(j,:) = y(j,:)/sum(y(j,:));
	end
	errorShade(100*xx,mean(y),sem(y),'Color',c(cc,:),'LineWidth',5);
	cc = cc+1;
end
xlabel('Flow Rate (mL/min)')
subplot(2,1,2), hold on
cc = 1;
for i = haz_data
	hist_this = data(i).PID;
	xx =  linspace(min(min(hist_this)),max(max(hist_this)),50);
	y = NaN(width(hist_this),50);
	for j = 1:width(hist_this)
		y(j,:) = hist(hist_this(:,j),xx);
		y(j,:) = y(j,:)/sum(y(j,:));
	end
	errorShade(xx,mean(y),sem(y),'Color',c(cc,:),'LineWidth',5);
	cc = cc+1;
end
xlabel('Stimulus(V)')
PrettyFig()

if being_published
	snapnow
	delete(gcf)
end

%%
% Clearly, there is some weird transformation between the flow rate through the odour and measured output. What if, instead of plotting the flow rates, we plot the fractional flow rate? 

figure('outerposition',[0 0 700 700],'PaperUnits','points','PaperSize',[700 700]); hold on
subplot(2,1,1), hold on
cc = 1;
for i = haz_data
	hist_this = data(i).MFC500;
	hist_this = hist_this*100./(hist_this*100 + 2000);
	xx =  linspace(min(min(hist_this)),max(max(hist_this)),50);
	y = NaN(width(hist_this),50);
	for j = 1:width(hist_this)
		y(j,:) = hist(hist_this(:,j),xx);
		y(j,:) = y(j,:)/sum(y(j,:));
	end
	errorShade(100*xx,mean(y),sem(y),'Color',c(cc,:),'LineWidth',5);
	cc = cc+1;
end
xlabel('Fractional Flow Rate (%)')
subplot(2,1,2), hold on
cc = 1;
for i = haz_data
	hist_this = data(i).PID;
	xx =  linspace(min(min(hist_this)),max(max(hist_this)),50);
	y = NaN(width(hist_this),50);
	for j = 1:width(hist_this)
		y(j,:) = hist(hist_this(:,j),xx);
		y(j,:) = y(j,:)/sum(y(j,:));
	end
	errorShade(xx,mean(y),sem(y),'Color',c(cc,:),'LineWidth',5);
	cc = cc+1;
end
xlabel('Stimulus(V)')
PrettyFig()

if being_published
	snapnow
	delete(gcf)
end

%%
% It's still weird, and there is still a bizarre non-linearity. Can we fit this with a  general purpose NLN model? 

fit_me = data(haz_data);
for i = 1:length(fit_me)
	fit_me(i).stimulus = mean2(fit_me(i).MFC500);
	fit_me(i).stimulus = (100*fit_me(i).stimulus)./(100*fit_me(i).stimulus+2000);
	fit_me(i).response = mean2(fit_me(i).PID);
end

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

t = toc;

%% 
% This document was built in: 
disp(strcat(oval(t,3),' seconds.'))

% tag the file as being published 
% add homebrew path
path1 = getenv('PATH');
path1 = [path1 ':/usr/local/bin'];
setenv('PATH', path1);

if being_published
	unix(strjoin({'tag -a published',which(mfilename)}));
end
