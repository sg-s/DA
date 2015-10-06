% makeFig9.m
% makes figure 9, which shows gain changes even when driven with light through reachr
% 
% created by Srinivas Gorur-Shandilya at 9:50 , 06 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


% add homebrew path
path1 = getenv('PATH');
if isempty(strfind(path1,[':/usr/local/bin']))
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


figure('outerposition',[0 0 1500 800],'PaperUnits','points','PaperSize',[1500 800]); hold on

axes_handles(1) = subplot(3,5,2:3);
axes_handles(2) = subplot(3,5,4);
axes_handles(3) = subplot(3,5,5);

axes_handles(4) = subplot(3,5,7:8);
axes_handles(5) = subplot(3,5,9);
axes_handles(6) = subplot(3,5,10); 

axes_handles(7) = subplot(3,5,12:13);
axes_handles(8) = subplot(3,5,14);
axes_handles(9) = subplot(3,5,15); 

for i = 1:length(axes_handles)
	hold(axes_handles(i),'on');
end 

%  #######  ########   #######  ########           #######  ########   #######  ########  
% ##     ## ##     ## ##     ## ##     ##         ##     ## ##     ## ##     ## ##     ## 
% ##     ## ##     ## ##     ## ##     ##         ##     ## ##     ## ##     ## ##     ## 
% ##     ## ##     ## ##     ## ########  ####### ##     ## ##     ## ##     ## ########  
% ##     ## ##     ## ##     ## ##   ##           ##     ## ##     ## ##     ## ##   ##   
% ##     ## ##     ## ##     ## ##    ##          ##     ## ##     ## ##     ## ##    ##  
%  #######  ########   #######  ##     ##          #######  ########   #######  ##     ## 

% load the MSG data
% load cached data
load('../data/MeanShiftedGaussians.mat')

% shorten paradigm names by throwing out 'MFC'
short_paradigm_names = paradigm_names;
for i = 1:length(paradigm_names)
	short_paradigm_names{i} = paradigm_names{i}(strfind(paradigm_names{i},'-')+1:end);
end
load('../data/MSG_per_neuron.mat','MSG_data')

c = parula(9); % colormap
ss = 50;
x = mean2([MSG_data(1,:).resp]);
y = mean2([MSG_data(5,:).resp]);
plot(axes_handles(1),MSG_data(1,1).time(1:ss:end),x(1:ss:end),'Color',c(1,:));
plot(axes_handles(1),MSG_data(1,1).time(1:ss:end),y(1:ss:end),'Color',c(8,:));
plot(axes_handles(2),[0 45],[0 45],'k--')
plot(axes_handles(2),x(1:ss:end),y(1:ss:end),'.','Color',c(8,:));
set(axes_handles(2),'XLim',[0 45],'YLim',[0 45])
set(axes_handles(1),'XLim',[35 55])

%  #######  ########   #######  ########          ##       ####  ######   ##     ## ######## 
% ##     ## ##     ## ##     ## ##     ##         ##        ##  ##    ##  ##     ##    ##    
% ##     ## ##     ## ##     ## ##     ##         ##        ##  ##        ##     ##    ##    
% ##     ## ##     ## ##     ## ########  ####### ##        ##  ##   #### #########    ##    
% ##     ## ##     ## ##     ## ##   ##           ##        ##  ##    ##  ##     ##    ##    
% ##     ## ##     ## ##     ## ##    ##          ##        ##  ##    ##  ##     ##    ##    
%  #######  ########   #######  ##     ##         ######## ####  ######   ##     ##    ##    


[PID, LFP, fA, paradigm, orn, fly, AllControlParadigms, paradigm_hashes, sequence] = consolidateData('/local-data/DA-paper/reachr/odour-flicker-light-background/',true);

% recompute paradigms sensibly 
light_V = NaN*paradigm;
for i = 1:length(paradigm)
	p = paradigm(i);
	if isempty(strfind(AllControlParadigms(p).Name,'LightFlicker')) & ~isempty(strfind(AllControlParadigms(p).Name,'Odour')) & ~isempty(strfind(AllControlParadigms(p).Name,'Light'))
		light_V(i) = str2double(AllControlParadigms(p).Name(strfind(AllControlParadigms(p).Name,'+')+1:strfind(AllControlParadigms(p).Name,'V')-1));
	end
end


x = mean2(fA(:,light_V == 0));
y = mean2(fA(:,light_V == 3.5));

time = 1e-3*(1:length(x));
plot(axes_handles(4),time,x,'Color',c(1,:))
plot(axes_handles(4),time,y,'Color',c(8,:))
set(axes_handles(4),'XLim',[35 55])

plot(axes_handles(5),[0 45],[0 45],'k--')
plot(axes_handles(5),x(20e3:ss:end),y(20e3:ss:end),'.','Color',c(8,:))
set(axes_handles(5),'XLim',[0 45],'YLim',[0 45])

% remove some junk
light_V(light_V == 4) = NaN;
light_V(light_V == .8) = NaN;
light_V(light_V == 1.1) = NaN;
light_V(light_V == 1.2) = NaN;

% compute gains for each light case, per trial
gain = NaN*paradigm;
for i = 1:length(paradigm)
	y = fA(:,i);
	try
		ff = fit(x(20e3:end),y(20e3:end),'Poly1');
		gain(i) = ff.p1;
	end
end

all_light_V = nonnans(unique(light_V));
all_gain = NaN*all_light_V;
all_gain_err = all_gain;
n = all_gain;

for i = 1:length(all_light_V)
	n(i) = sum(light_V == all_light_V(i));

	all_gain(i) = mean(gain(light_V == all_light_V(i)));
	all_gain_err(i) = std(gain(light_V == all_light_V(i)));
	all_gain_err(i) = all_gain_err(i)/sqrt(n(i));
end

axes(axes_handles(6))
errorbar(all_light_V,all_gain,all_gain_err,'k.')


% ##       ####  ######   ##     ## ########          #######  ########   #######  ########  
% ##        ##  ##    ##  ##     ##    ##            ##     ## ##     ## ##     ## ##     ## 
% ##        ##  ##        ##     ##    ##            ##     ## ##     ## ##     ## ##     ## 
% ##        ##  ##   #### #########    ##    ####### ##     ## ##     ## ##     ## ########  
% ##        ##  ##    ##  ##     ##    ##            ##     ## ##     ## ##     ## ##   ##   
% ##        ##  ##    ##  ##     ##    ##            ##     ## ##     ## ##     ## ##    ##  
% ######## ####  ######   ##     ##    ##             #######  ########   #######  ##     ## 


[PID, LFP, fA, paradigm, orn, fly, AllControlParadigms, paradigm_hashes, sequence] = consolidateData('/local-data/DA-paper/reachr/light-flicker-odour-background/',true);


% recompute paradigms sensibly 
odour_background = NaN*paradigm;
odour_V = NaN*paradigm;
for i = 1:length(paradigm)
	p = paradigm(i);
	if isempty(strfind(AllControlParadigms(p).Name,'OdourFlicker')) & ~isempty(strfind(AllControlParadigms(p).Name,'Odour')) & ~isempty(strfind(AllControlParadigms(p).Name,'LightFlicker'))
		odour_background(i) = mean(PID(20e3:end,i));

		odour_V(i) = str2double(AllControlParadigms(p).Name(strfind(AllControlParadigms(p).Name,'+')+1:strfind(AllControlParadigms(p).Name,'V')-1));
	end
end


x = mean2(fA(:,odour_V ==  0));
y = mean2(fA(:,odour_V == .5));

time = 1e-3*(1:length(x));
plot(axes_handles(7),time,x,'Color',c(1,:))
plot(axes_handles(7),time,y,'Color',c(8,:))
set(axes_handles(7),'XLim',[35 55])

plot(axes_handles(8),[0 45],[0 45],'k--')
plot(axes_handles(8),x(20e3:ss:end),y(20e3:ss:end),'.','Color',c(8,:))
set(axes_handles(8),'XLim',[0 45],'YLim',[0 45])

% compute gains for each light case, per trial
gain = NaN*paradigm;
for i = 1:length(paradigm)
	y = fA(:,i);
	try
		ff = fit(x(20e3:end),y(20e3:end),'Poly1');
		gain(i) = ff.p1;
	end
end

all_odour_V = nonnans(unique(odour_V));
all_gain = NaN*all_odour_V;
all_gain_err = all_gain;
n = all_gain;

for i = 1:length(all_odour_V)
	n(i) = sum(odour_V == all_odour_V(i));

	all_gain(i) = nanmean(gain(odour_V == all_odour_V(i)));
	all_gain_err(i) = nanstd(gain(odour_V == all_odour_V(i)));
	all_gain_err(i) = all_gain_err(i)/sqrt(n(i));
end

% convert control odour voltage to measured odour voltage
measured_odour = [];
for i = 1:length(all_odour_V)
	measured_odour(i) = mean(mean(PID(20e3:end,odour_V == all_odour_V(i))));
end

axes(axes_handles(9))
errorbar(measured_odour(n>4),all_gain(n>4),all_gain_err(n>4),'k.')

% cosmetics, labels, etc. 
ylabel(axes_handles(1),'Response (Hz)')

xlabel(axes_handles(2),'Response (Hz)')
ylabel(axes_handles(2),['Response to' char(10) 'stimulus + background (Hz)'])

ylabel(axes_handles(3),'Relative Gain')
xlabel(axes_handles(3),'Background odour (V)')

ylabel(axes_handles(4),'Response (Hz)')

xlabel(axes_handles(5),'Response (Hz)')
ylabel(axes_handles(5),['Response to' char(10) 'stimulus + background (Hz)'])



set(axes_handles(6),'XLim',[-.1 4],'YLim',[.2 1.2])
ylabel(axes_handles(6),'Relative Gain')
xlabel(axes_handles(6),'Background light (V)')

xlabel(axes_handles(7),'Time (s)')
ylabel(axes_handles(7),'Response (Hz)')

xlabel(axes_handles(8),'Response (Hz)')
ylabel(axes_handles(8),['Response to' char(10) 'stimulus + background (Hz)'])

set(axes_handles(9),'XLim',[-.01 .6],'YLim',[.2 1.2])
ylabel(axes_handles(9),'Relative Gain')
xlabel(axes_handles(9),'Background odour (V)')

prettyFig('fs=15;')

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
