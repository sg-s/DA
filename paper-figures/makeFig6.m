% makeFig6
% makes figure 6 of the paper, showing that gain control is broadly observed
% 
% created by Srinivas Gorur-Shandilya at 9:56 , 04 October 2015. Contact me at http://srinivas.gs/contact/
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

%      ########     ###             ##      ## #### ########  ######## ##       ##    ## 
%      ##     ##   ## ##            ##  ##  ##  ##  ##     ## ##       ##        ##  ##  
%      ##     ##  ##   ##           ##  ##  ##  ##  ##     ## ##       ##         ####   
%      ##     ## ##     ##          ##  ##  ##  ##  ##     ## ######   ##          ##    
%      ##     ## #########          ##  ##  ##  ##  ##     ## ##       ##          ##    
%      ##     ## ##     ##          ##  ##  ##  ##  ##     ## ##       ##          ##    
%      ########  ##     ##           ###  ###  #### ########  ######## ########    ## 

%       #######  ########   ######  ######## ########  ##     ## ######## ########  
%      ##     ## ##     ## ##    ## ##       ##     ## ##     ## ##       ##     ## 
%      ##     ## ##     ## ##       ##       ##     ## ##     ## ##       ##     ## 
%      ##     ## ########   ######  ######   ########  ##     ## ######   ##     ## 
%      ##     ## ##     ##       ## ##       ##   ##    ##   ##  ##       ##     ## 
%      ##     ## ##     ## ##    ## ##       ##    ##    ## ##   ##       ##     ## 
%       #######  ########   ######  ######## ##     ##    ###    ######## ########  


%% Figure 6: Gain Control is widely observed
% The previous result showed fast gain control using ethyl acetate, a volatile odor, in ab3A, an antennal ORN. To determine if fast gain control is widely observed, and not a peculiarity of the odor-receptor combination used, I reanalyzed a published data set of ORN responses to flickering odor stimuli (from Martelli et al. J. Neuro. 2013). The data set consists of six different odors and two ORN types: one in the antenna and one in the maxillary palp. My analysis of fast gain control is self-consistent across experimental replicates (A-C), and shows that fast gain control is observed for different odors (D-F) and in two different neurons (G-I).  In each row, the first column shows the linear filter extracted from the stimulus and the response. The second column shows the principal components fit to the highest (lowest) 1/3 of stimulus filtered over some history length in red (green). The third column shows how relative gain varies at times when the filtered stimulus in in the top 1/3 (green) or bottom 1/3 (red), for various history lengths. Only significant (p<0.01) points are shown. In every dataset, significant gain control is observed on a sub-second time scale, and the gain control always acts to amplify responses to recently weak stimuli and suppress responses to recently large stimuli (red dots always below green dots in third column). 

load('../data/CMData_Gain.mat')
load('../data/CM_Data_filters.mat')
combined_data_file = ('/local-data/DA-paper/carlotta-martelli/flickering-stim/data.mat');
load(combined_data_file)
filtertime = -200:700;
filtertime = filtertime*1e-3;

hl_min = .2;
hl_max = 10;
history_lengths = [logspace(log10(hl_min),log10(.5),15) logspace(log10(.5),log10(10),15)];
history_lengths = unique(history_lengths);


s = 860;
figure('outerposition',[0 0 s s],'PaperUnits','points','PaperSize',[s s]); hold on, clear s
clear axes_handles
for i = 1:9
	axes_handles(i) = subplot(3,3,i); hold on
end

% ########  ######## ########  ##       ####  ######     ###    ######## ########  ######  
% ##     ## ##       ##     ## ##        ##  ##    ##   ## ##      ##    ##       ##    ## 
% ##     ## ##       ##     ## ##        ##  ##        ##   ##     ##    ##       ##       
% ########  ######   ########  ##        ##  ##       ##     ##    ##    ######    ######  
% ##   ##   ##       ##        ##        ##  ##       #########    ##    ##             ## 
% ##    ##  ##       ##        ##        ##  ##    ## ##     ##    ##    ##       ##    ## 
% ##     ## ######## ##        ######## ####  ######  ##     ##    ##    ########  ######  

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


% first row: experimental replicates
do_these = [2 7 9 13 15 16];

for i = do_these
	K = allfilters(i).K;
	for j = 1:width(K)
		K(:,j) = K(:,j)/max(K(:,j));
	end
	plot(axes_handles(1),filtertime,mean(K,2))
end
clear ph
ph(3:4) = axes_handles(2:3);


for i = do_these
	time = 1e-3*(1:length(mean(data(i).PID,2)));
	stimulus = mean(data(i).PID,2);
	prediction = mean(data(i).LinearFit,2);
	response = mean(data(i).fA,2);

	% throw out first 5 seconds
	time = time(5e3:end);
	stimulus = stimulus(5e3:end);
	response = response(5e3:end);
	prediction = prediction(5e3:end);

	% remove trend in stimulus
	temp = fit(time(:),stimulus(:),'poly2');
	stimulus = stimulus - temp(time) + mean(stimulus);

	% fix the gain to be exactly 1
	x = prediction;
	y = response;
	rm_this = isnan(x) | isnan(y) ;
	x(rm_this) = [];
	y(rm_this) = [];
	temp = fit(x,y,'poly1');
	prediction = prediction*temp.p1;

	% ignore very low responses
	y(y<5) = NaN;

	gainAnalysisWrapper('time',time,'response',response,'stimulus',stimulus,'prediction',prediction,'history_lengths',history_lengths,'ph',ph,'engine',@gainAnalysis,'use_cache',true,'example_history_length',history_lengths(15));
end


% add a minimal legend
h=get(ph(3),'Children');
legend(h(1:2),{'High Stim.','Low Stim.'},'Location','northwest')

%    ########  #### ######## ######## ######## ########  ######## ##    ## ######## 
%    ##     ##  ##  ##       ##       ##       ##     ## ##       ###   ##    ##    
%    ##     ##  ##  ##       ##       ##       ##     ## ##       ####  ##    ##    
%    ##     ##  ##  ######   ######   ######   ########  ######   ## ## ##    ##    
%    ##     ##  ##  ##       ##       ##       ##   ##   ##       ##  ####    ##    
%    ##     ##  ##  ##       ##       ##       ##    ##  ##       ##   ###    ##    
%    ########  #### ##       ##       ######## ##     ## ######## ##    ##    ##    
   
%    ##    ## ######## ##     ## ########   #######  ##    ##  ######  
%    ###   ## ##       ##     ## ##     ## ##     ## ###   ## ##    ## 
%    ####  ## ##       ##     ## ##     ## ##     ## ####  ## ##       
%    ## ## ## ######   ##     ## ########  ##     ## ## ## ##  ######  
%    ##  #### ##       ##     ## ##   ##   ##     ## ##  ####       ## 
%    ##   ### ##       ##     ## ##    ##  ##     ## ##   ### ##    ## 
%    ##    ## ########  #######  ##     ##  #######  ##    ##  ######  

% add new ab2 data
p = '/local-data/DA-paper/large-variance-flicker/ab2';
[PID, ~, fA, paradigm, orn, fly] = consolidateData(p,true);


do_these = [11 17];
l = [];
for i = do_these
	K = allfilters(i).K;
	for j = 1:width(K)
		K(:,j) = K(:,j)/max(K(:,j));
	end
	l=[l plot(axes_handles(7),filtertime,mean(K,2))];
end
legend(l,{'pb1A','ab3A'})


clear ph
ph(3:4) = axes_handles(8:9);
for i = do_these
	time = 1e-3*(1:length(mean(data(i).PID,2)));
	stimulus = mean(data(i).PID,2);
	prediction = mean(data(i).LinearFit,2);
	response = mean(data(i).fA,2);

	% throw out first 5 seconds
	time = time(5e3:end);
	stimulus = stimulus(5e3:end);
	response = response(5e3:end);
	prediction = prediction(5e3:end);

	% remove trend in stimulus
	temp = fit(time(:),stimulus(:),'poly2');
	stimulus = stimulus - temp(time) + mean(stimulus);

	% fix the gain to be exactly 1
	x = prediction;
	y = response;
	rm_this = isnan(x) | isnan(y) ;
	x(rm_this) = [];
	y(rm_this) = [];
	temp = fit(x,y,'poly1');
	prediction = prediction*temp.p1;

	% ignore very low responses
	y(y<5) = NaN;

	gainAnalysisWrapper('time',time,'response',response,'stimulus',stimulus,'prediction',prediction,'history_lengths',history_lengths,'ph',ph,'engine',@gainAnalysis,'use_cache',true,'example_history_length',history_lengths(15)');
end



% add a minimal legend
h=get(ph(3),'Children');
legend(h(1:2),{'High Stim.','Low Stim.'},'Location','northwest')


%     ########  #### ######## ######## ######## ########  ######## ##    ## ######## 
%     ##     ##  ##  ##       ##       ##       ##     ## ##       ###   ##    ##    
%     ##     ##  ##  ##       ##       ##       ##     ## ##       ####  ##    ##    
%     ##     ##  ##  ######   ######   ######   ########  ######   ## ## ##    ##    
%     ##     ##  ##  ##       ##       ##       ##   ##   ##       ##  ####    ##    
%     ##     ##  ##  ##       ##       ##       ##    ##  ##       ##   ###    ##    
%     ########  #### ##       ##       ######## ##     ## ######## ##    ##    ##    
    
%      #######  ########   #######  ##     ## ########   ######  
%     ##     ## ##     ## ##     ## ##     ## ##     ## ##    ## 
%     ##     ## ##     ## ##     ## ##     ## ##     ## ##       
%     ##     ## ##     ## ##     ## ##     ## ########   ######  
%     ##     ## ##     ## ##     ## ##     ## ##   ##         ## 
%     ##     ## ##     ## ##     ## ##     ## ##    ##  ##    ## 
%      #######  ########   #######   #######  ##     ##  ######  


do_these = [7 8 10 14 17 18];
odours = {'1but','1o3ol','dsucc','2ac','2but','5ol'};
l = [];
for i = do_these
	K = allfilters(i).K;
	for j = 1:width(K)
		K(:,j) = K(:,j)/max(K(:,j));
	end
	l=[l plot(axes_handles(4),filtertime,mean(K,2))];
end
legend(l,odours)

clear ph
ph(3:4) = axes_handles(5:6);
for i = do_these
	time = 1e-3*(1:length(mean(data(i).PID,2)));
	stimulus = mean(data(i).PID,2);
	prediction = mean(data(i).LinearFit,2);
	response = mean(data(i).fA,2);

	% throw out first 5 seconds
	time = time(5e3:end);
	stimulus = stimulus(5e3:end);
	response = response(5e3:end);
	prediction = prediction(5e3:end);

	% remove trend in stimulus
	temp = fit(time(:),stimulus(:),'poly2');
	stimulus = stimulus - temp(time) + mean(stimulus);

	% fix the gain to be exactly 1
	x = prediction;
	y = response;
	rm_this = isnan(x) | isnan(y) ;
	x(rm_this) = [];
	y(rm_this) = [];
	temp = fit(x,y,'poly1');
	prediction = prediction*temp.p1;

	gainAnalysisWrapper('time',time,'response',response,'stimulus',stimulus,'prediction',prediction,'history_lengths',history_lengths,'ph',ph,'engine',@gainAnalysis,'use_cache',true,'example_history_length',history_lengths(15));
end


% add a minimal legend
h=get(ph(3),'Children');
legend(h(1:2),{'High Stim.','Low Stim.'},'Location','northwest')

% clean up -- remove the scatter points
for j = [2 5 8]
	h=get(axes_handles(j),'Children');
	rm_this = [];
	for i = 1:length(h)
		if strcmp(get(h(i),'Marker'),'.')
			rm_this = [rm_this i];
		end
	end
	delete(h(rm_this))
end

% remove the line indicating the example history plot
for j = [3 6 9]
	h=get(axes_handles(j),'Children');
	rm_this = [];
	for i = 1:length(h)
		try
			if  strcmp(get(h(i),'LineStyle'),'-.')
				rm_this = [rm_this i];
			end
		catch
		end
	end
	delete(h(rm_this))

	% remove the lines where the data isn't significant
	h=get(axes_handles(j),'Children');
	rm_this = [];
	for i = 1:length(h)
		try
			if  strcmp(get(h(i),'Type'),'line')
				rm_this = [rm_this i];
			end
		catch
		end
	end
	delete(h(rm_this))

	% make all the scatter plots smaller
	h=get(axes_handles(j),'Children');
	for i = 1:length(h)
		 set(h(i),'SizeData',256)
	end

	set(axes_handles(j),'XLim',[.1 10],'YLim',[.6 1.45])
end


% cosmetics
xlabel(axes_handles(1),'Lag (s)')
xlabel(axes_handles(4),'Lag (s)')
title(axes_handles(1),'Filters')
title(axes_handles(4),'Filters')
title(axes_handles(7),'Filters')
ylabel(axes_handles(3),'Relative Gain')
ylabel(axes_handles(6),'Relative Gain')
ylabel(axes_handles(9),'Relative Gain')
title(axes_handles(2),'')
title(axes_handles(5),'')
title(axes_handles(8),'')



prettyFig('plw=1.5;','lw=1.5;','fs=14;')

ylabel(axes_handles(1),'Exp. Replicates','FontSize',20)
ylabel(axes_handles(7),'Diff. ORNs','FontSize',20)
ylabel(axes_handles(4),'Diff. odors','FontSize',20)

if being_published
	snapnow
	delete(gcf)
end

%% Pulse Gain Analysis
% In this section we analyse the gain in a differnet way. We compute the gain of the ORN at every valve onset, and try to find a projection of the stimulus that accounts for the variation in this gain. 

history_lengths = logspace(log10(.1),log10(10),40);
gain_K1_slope = NaN(length(history_lengths),length(data));
gain_K1_rho = NaN(length(history_lengths),length(data));

for i = 1:length(data)


	valve = data(i).Valve';
	inst_gain = findGainWhenValveOpens(valve,stim,resp);
	for j = 1:length(history_lengths)

		% first do the simple box filter
		filtered_stim = floor(history_lengths(j)*1e3);
		filtered_stim = filter(ones(filtered_stim,1),filtered_stim,stim);

		y = inst_gain;
		y(1:10e3) = NaN;

		rm_this = isnan(filtered_stim) | isnan(y);
		filtered_stim(rm_this) = [];
		y(rm_this) = [];

		ff = fit(filtered_stim(:),y(:),'power1','StartPoint',[0 0]);
		gain_K1_rho(j,i) = spear(filtered_stim(1:10:end),y(1:10:end));
		gain_K1_slope(j,i) = ff.b;
	end
end




%% Instantaneous gain analysis
% In this section we analyze the data in a different way: we first compute the instantaneous gain for all the data, and try to find projections of the stimulus that maximally predict the instantaneous gain. 

% compute the instantenous gain for each case

gain_K1_slope = NaN(length(history_lengths),length(data));
gain_K1_rho = NaN(length(history_lengths),length(data));

gain_K2_slope = NaN(length(history_lengths),length(data));
gain_K2_rho = NaN(length(history_lengths),length(data));

gain_K3_slope = NaN(length(history_lengths),length(data));
gain_K3_rho = NaN(length(history_lengths),length(data));


for i = 1:length(data)

	t = 1e-3*(1:length(data(i).PID));
	filtertime = 1e-3*(1:length(best_reg_filters(i).K)) - .2;

	stim = data(i).PID;
	resp = data(i).fA;

	stim = mean(stim,2);
	resp = mean(resp,2);

	pred = mean(stim) + convolve(t,stim,best_reg_filters(i).K,filtertime);

	[~,inst_gain] = makeFig6G(stim,resp,pred,500);

	for j = 1:length(history_lengths)
		stim = mean(data(i).PID,2);
		resp = mean(data(i).fA,2);

		y = inst_gain;
		y(1:10e3) = -1;

		% first do the simple box filter
		temp = floor(history_lengths(j)*1e3);
		temp = filter(ones(temp,1),temp,stim);

		rm_this = (isnan(temp) | isnan(inst_gain) | y < 0);
		temp(rm_this) = [];
		y(rm_this) = [];

		ff = fit(temp(:),y(:),'power1','StartPoint',[0 0]);
		temp = ff(temp);
		gain_K1_rho(j,i) = spear(temp(1:10:end),y(1:10:end));
		gain_K1_slope(j,i) = ff.b;
	end
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
