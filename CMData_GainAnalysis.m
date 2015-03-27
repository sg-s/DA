% CMData_GainAnalysis.m
% 
% created by Srinivas Gorur-Shandilya at 2:31 , 26 March 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

%%
% In this document we look at carlott's data from a gain analysis perspective. 

%% Valve switching time affects gain properties 
% The first interesting thing we see in this data is that the switching time of the valve affects the overall picture we get from gain analysis. In the following figure, we categorise all the data by the switching time of the valve. 

%%
% As can be seen from the figure, when the valve switching time is < 100ms, we don't consistently see the green above the red anywhere. Only when the valve switching time is > 100ms, do we see a region where the green lines are above the red. 


calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end
tic

load('CMData_Gain.mat')
load('CM_Data_filters.mat')
combined_data_file = ('/local-data/DA-paper/carlotta-martelli/flickering-stim/data.mat');
load(combined_data_file)
filtertime = -200:700;
filtertime = filtertime*1e-3;

% add a new field with correlation times
for i = 1:length(data)
	if any(strfind(data(i).original_name,'30ms'))
		data(i).corr_time = 30;
	end
	if any(strfind(data(i).original_name,'50ms'))
		data(i).corr_time = 50;
	end
	if any(strfind(data(i).original_name,'100ms'))
		data(i).corr_time = 100;
	end
end


figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
for i = 1:3
	ax(i) = subplot(1,3,i); hold on
	set(gca,'XScale','log')
end

title(ax(1),'30ms correlation')
title(ax(2),'50ms correlation')
title(ax(3),'100ms correlation')

history_lengths = (3*floor(1000*logspace(-1,1,30)/3))/1e3;
p = NaN(2,length(history_lengths));

cross_over_time = NaN*history_lengths;
filter_peak_loc = NaN*history_lengths;
filter_width = NaN*history_lengths;
stimulus_autocorr = NaN*history_lengths;
response_autocorr = NaN*history_lengths;

all_low = NaN(length(data),30);
all_high = NaN(length(data),30);

for i = 1:length(data)
	% figure out where to plot
	if any(strfind(data(i).original_name,'30ms'))
		plot_here = ax(1);
	end
	if any(strfind(data(i).original_name,'50ms'))
		plot_here = ax(2);
	end
	if any(strfind(data(i).original_name,'100ms'))
		plot_here = ax(3);
	end

	ph = [];
	ph(4) = plot_here;

	IGA_data.time = 1e-3*(1:length(mean2(data(i).PID)));
	IGA_data.stimulus = mean2(data(i).PID);
	IGA_data.prediction = mean2(data(i).LinearFit);
	IGA_data.response = mean2(data(i).fA);

	% throw out first 5 seconds
	IGA_data.time = IGA_data.time(5e3:end);
	IGA_data.stimulus = IGA_data.stimulus(5e3:end);
	IGA_data.response = IGA_data.response(5e3:end);
	IGA_data.prediction = IGA_data.prediction(5e3:end);

	IGA_data.frac = .33;

	% remove trend in stimulus
	temp = fit(IGA_data.time(:),IGA_data.stimulus(:),'poly2');
	IGA_data.stimulus = IGA_data.stimulus - temp(IGA_data.time) + mean(IGA_data.stimulus);

	% add the name
	IGA_data.name = data(i).original_name;

	% fix the gain to be exactly 1
	x = IGA_data.prediction;
	y = IGA_data.response;
	rm_this = isnan(x) | isnan(y);
	x(rm_this) = [];
	y(rm_this) = [];
	temp = fit(x,y,'poly1');
	IGA_data.prediction = IGA_data.prediction*temp.p1;



	[~,low_slopes,high_slopes]=GainAnalysis4(IGA_data,history_lengths,[],ph,p);

	all_low(i,:) = low_slopes;
	all_high(i,:) = high_slopes;

	if any(strfind(data(i).original_name,'100ms'))
		temp = find(high_slopes - low_slopes < 0,1,'first');
		if isempty(temp)
			temp = NaN;
		else
			cross_over_time(i) = history_lengths(temp);
			K = mean2(allfilters(i).K);
			[~,loc]=max(K);
			filter_peak_loc(i) = loc-200;
			filter_width(i) = find(abs(K(200:end))>3*std(K(1:200)),1,'last');
			[~,~,~,x4]=FindCorrelationTime(IGA_data.stimulus);
			stimulus_autocorr(i) = x4;
			[~,~,~,x4]=FindCorrelationTime(IGA_data.response);
			response_autocorr(i) = x4;


		end
	end
end


PrettyFig('plw=1.5;','lw=1.5;','fs=14;')

if being_published
	snapnow
	delete(gcf)
end


%% What determines cross over time? 
% What's also interesting is that when the valve switching time is > 100ms, the green lines and the red lines switch over at some point. What determines this switch over time? 

%%
% The following figure shows all the 100ms correlated data, grouped by odour type, with dots to indicate significant differences:

figure('outerposition',[0 0 1500 1000],'PaperUnits','points','PaperSize',[1500 1000]); hold on
for j = 1:6
	ax(j) = subplot(2,3,j); hold on
	set(gca,'XScale','log')
end

odours = {'1but','2but','2ac','5ol','1o3ol','d2succ'};

for i = 1:length(data)
	% figure out where to plot
	plot_here = [];
	for j = 1:length(odours)
		if any(strfind(data(i).original_name,'100ms')) && any(strfind(data(i).original_name,odours{j}))
			plot_here = ax(j);
		end
	end
	
	if ~isempty(plot_here)
		clear x
		x.time = 1e-3*(1:length(mean2(data(i).PID)));
		x.stimulus = mean2(data(i).PID);
		x.prediction = mean2(data(i).LinearFit);
		x.response = mean2(data(i).fA);

		% throw out first 5 seconds
		x.time = x.time(5e3:end);
		x.stimulus = x.stimulus(5e3:end);
		x.response = x.response(5e3:end);
		x.prediction = x.prediction(5e3:end);

		% remove trend in stimulus
		temp = fit(x.time(:),x.stimulus(:),'poly2');
		x.stimulus = x.stimulus - temp(x.time) + mean(x.stimulus);

		% fix the gain to be exactly 1 -- this is trivial -- the gain is not one because we average over many trials. 
		a = x.prediction;
		b = x.response;
		rm_this = isnan(a) | isnan(b);
		a(rm_this) = [];
		b(rm_this) = [];
		temp = fit(a,b,'poly1');
		x.prediction = x.prediction*temp.p1;


		x.frac = .33;

		clear ph
		ph(4) = plot_here;
		[p,~,~,~,~,history_lengths]=GainAnalysisWrapper2('response',x.response,'prediction',x.prediction,'stimulus',x.stimulus,'time',x.time,'ph',ph,'history_lengths',history_lengths);
	end
end

for j = 1:6
	title(ax(j),odours{j})
end

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%%
% That's interesting. We will look in detail at one of the cases where this switch over happens:
% 
% <</code/da/CMData_fig1.PNG>>
%

%%
% It looks like this is entirely because the neuron stops firing. What if we ignore all the times when the neuron is silent (f<5Hz)? 
% 
% <</code/da/CMData_fig2.PNG>>
%

%%
% Armed with this information, we can go back to investigating why we don't see the standard picture when the stimulus is correlated on 30 or 50ms.

%%
% I have no idea, but for 30ms and 50ms correlated stimuli, gain analysis seems to suggest that the fast gain control is inverted (gain at times when stimulus is locally high is HIGH), and this is not due to trivial reasons like neuron silencing, etc. It also appears to be odour-independent. 
% 
% <</code/da/CMData_fig3.PNG>>
%



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


