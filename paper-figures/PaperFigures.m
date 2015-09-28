% Paper Figures
% makes all the figures for the paper
% 
% created by Srinivas Gorur-Shandilya at 12:57 , 21 January 2015. Contact me at http://srinivas.gs/contact/
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
tic

% this determines which figures to do. 
fig1 = true; 	% how we determine gain
fig2 = true;	% weber-like gain control
fig_supp1 = false;
fig3 = false;	% speed-gain tradeoff
fig4 = false;	% natualistic stimuli + fast gain control
fig5 = false;	% fast gain control widely observed
fig6 = false; 	% switching experiment
fig7 = false;	% LFP experiments
fig8 = false; 	% models to explain this
fig9 = false;


end

%   ######## ####  ######   ##     ## ########  ########    ######## 
%   ##        ##  ##    ##  ##     ## ##     ## ##          ##       
%   ##        ##  ##        ##     ## ##     ## ##          ##       
%   ######    ##  ##   #### ##     ## ########  ######      #######  
%   ##        ##  ##    ##  ##     ## ##   ##   ##                ## 
%   ##        ##  ##    ##  ##     ## ##    ##  ##          ##    ## 
%   ##       ####  ######    #######  ##     ## ########     ######  



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



if fig5

%% Figure 4: Gain Control is widely observed
% The previous result showed fast gain control using ethyl acetate, a volatile odor, in ab3A, an antennal ORN. To determine if fast gain control is widely observed, and not a peculiarity of the odor-receptor combination used, I reanalyzed a published data set of ORN responses to flickering odor stimuli (from Martelli et al. J. Neuro. 2013). The data set consists of six different odors and two ORN types: one in the antenna and one in the maxillary palp. My analysis of fast gain control is self-consistent across experimental replicates (A-C), and shows that fast gain control is observed for different odors (D-F) and in two different neurons (G-I).  In each row, the first column shows the linear filter extracted from the stimulus and the response. The second column shows the principal components fit to the highest (lowest) 1/3 of stimulus filtered over some history length in red (green). The third column shows how relative gain varies at times when the filtered stimulus in in the top 1/3 (green) or bottom 1/3 (red), for various history lengths. Only significant (p<0.01) points are shown. In every dataset, significant gain control is observed on a sub-second time scale, and the gain control always acts to amplify responses to recently weak stimuli and suppress responses to recently large stimuli (red dots always below green dots in third column). 

clearvars -except being_published fig*
load('CMData_Gain.mat')
load('CM_Data_filters.mat')
combined_data_file = ('/local-data/DA-paper/carlotta-martelli/flickering-stim/data.mat');
load(combined_data_file)
filtertime = -200:700;
filtertime = filtertime*1e-3;
history_lengths = logspace(log10(.200),1,30);


s = 860;
figure('outerposition',[0 0 s s],'PaperUnits','points','PaperSize',[s s]); hold on, clear s
axes_handles = [];
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
	plot(axes_handles(1),filtertime,mean2(K))
end
clear ph
ph(3:4) = axes_handles(2:3);


for i = do_these
	time = 1e-3*(1:length(mean2(data(i).PID)));
	stimulus = mean2(data(i).PID);
	prediction = mean2(data(i).LinearFit);
	response = mean2(data(i).fA);

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

	GainAnalysisWrapper2('time',time,'response',response,'stimulus',stimulus,'prediction',prediction,'history_lengths',history_lengths,'ph',ph,'engine',@GainAnalysis5,'use_cache',true);
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


do_these = [11 17];
l = [];
for i = do_these
	K = allfilters(i).K;
	for j = 1:width(K)
		K(:,j) = K(:,j)/max(K(:,j));
	end
	l=[l plot(axes_handles(7),filtertime,mean2(K))];
end
legend(l,{'pb1A','ab3A'})


clear ph
ph(3:4) = axes_handles(8:9);
for i = do_these
	time = 1e-3*(1:length(mean2(data(i).PID)));
	stimulus = mean2(data(i).PID);
	prediction = mean2(data(i).LinearFit);
	response = mean2(data(i).fA);

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

	GainAnalysisWrapper2('time',time,'response',response,'stimulus',stimulus,'prediction',prediction,'history_lengths',history_lengths,'ph',ph,'engine',@GainAnalysis5,'use_cache',true);
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
	l=[l plot(axes_handles(4),filtertime,mean2(K))];
end
legend(l,odours)

clear ph
ph(3:4) = axes_handles(5:6);
for i = do_these
	time = 1e-3*(1:length(mean2(data(i).PID)));
	stimulus = mean2(data(i).PID);
	prediction = mean2(data(i).LinearFit);
	response = mean2(data(i).fA);

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

	GainAnalysisWrapper2('time',time,'response',response,'stimulus',stimulus,'prediction',prediction,'history_lengths',history_lengths,'ph',ph,'engine',@GainAnalysis5,'use_cache',true);
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

end

%          ######## ####  ######   ##     ## ########  ########     #######  
%          ##        ##  ##    ##  ##     ## ##     ## ##          ##     ## 
%          ##        ##  ##        ##     ## ##     ## ##          ##        
%          ######    ##  ##   #### ##     ## ########  ######      ########  
%          ##        ##  ##    ##  ##     ## ##   ##   ##          ##     ## 
%          ##        ##  ##    ##  ##     ## ##    ##  ##          ##     ## 
%          ##       ####  ######    #######  ##     ## ########     #######  



%% Figure 6: ORNs change gain rapidly to naturalistic stimuli. 
% Odorant stimuli in nature are thought to comprise of brief pulses of odorant, occurring in clumps of whiffs, with very broady distributed intensities. We delivered a reproducible, well-controlled odorant stimulus that mimicked these statistics to ab3A neurons (A). This stimulus was very broadly distributed (B), and elicited responses from the neuron for every whiff. By comparing the response of the neuron to a linear model, we we estimated the gain of the neuron at each whiff, and plotted the response trajectories colour-coded by the mean stimulus in the preceding 500ms (C). Neuron gain scales as an inverse power law with the mean stimulus in the preceding 500ms, changing from whiff to whiff (D).   

clearvars -except being_published fig*

if fig6



end


%             ######## ####  ######   ##     ## ########  ########          ######## 
%             ##        ##  ##    ##  ##     ## ##     ## ##                ##    ## 
%             ##        ##  ##        ##     ## ##     ## ##                    ##   
%             ######    ##  ##   #### ##     ## ########  ######               ##    
%             ##        ##  ##    ##  ##     ## ##   ##   ##                  ##     
%             ##        ##  ##    ##  ##     ## ##    ##  ##                  ##     
%             ##       ####  ######    #######  ##     ## ########            ##     

%% Figure 7: Timescales of adaptation to step changes in mean and variance
% Stimulus with alternating periods of low mean and high mean, keeping the variance the same (A). Neuron responses were averaged over all periods, and plotted as time since the switch from low to high mean (B). We investigated the kinetics and gain of the neuron at different points relative to the switch (coloured bars in B). 
% 
% <</Users/sigbhu/code/da/images/fig7.png>>
%

%% Figure 8: Local Field Potential also shows evidence of fast and slow gain control.  
% 
% <</Users/sigbhu/code/da/images/fig8.png>>
%



%      ##     ##  #######  ########  ######## ##        ######  
%      ###   ### ##     ## ##     ## ##       ##       ##    ## 
%      #### #### ##     ## ##     ## ##       ##       ##       
%      ## ### ## ##     ## ##     ## ######   ##        ######  
%      ##     ## ##     ## ##     ## ##       ##             ## 
%      ##     ## ##     ## ##     ## ##       ##       ##    ## 
%      ##     ##  #######  ########  ######## ########  ######  

if fig9


%% Figure 9: Models to explain observed phenomena 
%

figure('outerposition',[0 0 1400 600],'PaperUnits','points','PaperSize',[1400 600]); hold on
clear axes_handles
for i = 1:10
	axes_handles(i) = subplot(2,5,i); hold on
end

% first column shows fit quality
for i = 1:5:10
	plot(axes_handles(i),[0 5],[0 5],'k--')
	xlabel(axes_handles(i),'(P_{S}/P_{N})^{1/2}','interpreter','tex')
	ylabel(axes_handles(i),'(P_{S}/P_{R})^{1/2}','interpreter','tex')
end

% Label the columns
title(axes_handles(1),'Fit Quality')
title(axes_handles(2),'Weber Scaling')
title(axes_handles(3),'Response Speedup')
title(axes_handles(4),'Pulse Speedup')
title(axes_handles(5),'Fast Gain Control')



%       ##       ##    ##    ##     ##  #######  ########  ######## ##       
%       ##       ###   ##    ###   ### ##     ## ##     ## ##       ##       
%       ##       ####  ##    #### #### ##     ## ##     ## ##       ##       
%       ##       ## ## ##    ## ### ## ##     ## ##     ## ######   ##       
%       ##       ##  ####    ##     ## ##     ## ##     ## ##       ##       
%       ##       ##   ###    ##     ## ##     ## ##     ## ##       ##       
%       ######## ##    ##    ##     ##  #######  ########  ######## ######## 

load('LVF_data.mat')
tA = 1e-3*(1:length(data(1).stim));
% extract filters for each neuron
for i = 1:length(data)
	[this_K, ~, filtertime_full] = FindBestFilter(mean2(data(i).stim),mean2(data(i).resp),[],'regmax=1;','regmin=1;','filter_length=1999;','offset=500;');
	filtertime_full = filtertime_full*1e-3;
	filtertime = 1e-3*(-200:900);
	this_K = interp1(filtertime_full,this_K,filtertime);
	data(i).K = this_K;
end

% make linear predictions and extract non-linearities 
for i = 1:length(data)
	fp = convolve(tA,mean2(data(i).stim),data(i).K,filtertime);
	R = mean2(data(i).resp);
	temp =fit(fp(~(isnan(fp) | isnan(R))),R(~(isnan(fp) | isnan(R))),'poly1');
	fp = fp*temp.p1;
	fp = fp+temp.p2;
	data(i).LinearFit = fp;
	data(i).LNFit = hill(p_LN(i),data(i).LinearFit);

	% show fit quality of DA model for Large Variance Flicker
	[qx, qy] = GeffenMeister(data(i).resp,data(i).LNFit);
	plot(axes_handles(1),qx,qy,'ko')
end

% show gain analysis -- da model
ph = []; ph(3) = axes_handles(5);
history_lengths = 0.4890;
time = 1e-3*(1:length(data(1).LinearFit));
p=GainAnalysisWrapper2('response',mean2([data.resp]),'prediction',mean2([data.LNFit]),'stimulus',mean2([data.stim]),'time',time,'ph',ph,'history_lengths',history_lengths,'example_history_length',history_lengths(1),'engine',@GainAnalysis5);


% show the p-value
axes(axes_handles(5))
if p(1) == 0
	text(10,60,'p<0.01')
else
	text(10,60,strkat('p = ',oval(p(1))))
end

%    ##       ##    ##    ########  ##     ## ##        ######  ########  ######  
%    ##       ###   ##    ##     ## ##     ## ##       ##    ## ##       ##    ## 
%    ##       ####  ##    ##     ## ##     ## ##       ##       ##       ##       
%    ##       ## ## ##    ########  ##     ## ##        ######  ######    ######  
%    ##       ##  ####    ##        ##     ## ##             ## ##             ## 
%    ##       ##   ###    ##        ##     ## ##       ##    ## ##       ##    ## 
%    ######## ##    ##    ##         #######  ########  ######  ########  ######  



load('PulseData.mat')

% average filters
K = mean2(vertcat(data.K));
K = K/max(K);

% create the data to be fit
% clear d
% time = 1e-3*(1:length(PulseData(1).stim(:,1)));
% for i = 1:length(PulseData)
% 	d(i).stimulus = convolve(time,mean2(PulseData(i).stim),K,filtertime);
% 	d(i).response = mean2(PulseData(i).resp);
% end

clear p
p.A =  225.1576;
p.k =  382.2206;
p.n =  0.3952;

% this stores quality metrics
qx = NaN(length(PulseData),1);
qy = NaN(length(PulseData),1);

% make predictions everywhere
for i = 1:length(PulseData)
	PulseData(i).pred = NaN*PulseData(i).stim;
	for j = 1:width(PulseData(i).stim)
		PulseData(i).pred(:,j) = convolve(time,PulseData(i).stim(:,j),K,filtertime);
		PulseData(i).pred(:,j) = hill(PulseData(i).pred(:,j),p);
	end
	[qx(i), qy(i)] = GeffenMeister(PulseData(i).resp,PulseData(i).pred);
end


% show this on the plot
plot(axes_handles(1),qx,qy,'k+')

% show the timing data
stim_half_time = NaN(1e4,1); % time it takes to go to half max
resp_half_time = NaN(1e4,1); % time it takes to go to half max
foreground_stim = NaN(1e4,1);
c = 1;
a = 1050;
z = 1500; % nominal stimulus start and stop
for i = 1:length(PulseData)
	for j = 1:width(PulseData(i).stim)
		this_stim = PulseData(i).stim(:,j);
		foreground_stim(c) = mean(this_stim(a:z));
		this_stim = this_stim/max(this_stim(a:z));
		this_resp = PulseData(i).resp(:,j);
		this_resp = this_resp -  mean(this_resp(1:a));
		this_resp = this_resp/max(this_resp(a:z));
		stim_half_time(c) = max([find(this_stim(a:z)>.5,1,'first') NaN]);
		resp_half_time(c) = max([find(this_resp(a:z)>.5,1,'first') NaN]);
		c = c+1;
	end
end
stim_half_time(c:end) = [];
resp_half_time(c:end) = [];
resp_time_data = resp_half_time - stim_half_time;
foreground_stim(c:end) = [];

% clean up
rm_this = foreground_stim < 1e-2 | resp_time_data < 0 | resp_time_data > 100;
foreground_stim(rm_this) = [];
resp_time_data(rm_this) = [];

plot(axes_handles(4),foreground_stim,resp_time_data,'k+')
set(axes_handles(4),'XScale','log','XLim',[1e-2 20],'YLim',[0 100])
xlabel(axes_handles(4),'Mean Stimulus (V)')
ylabel(axes_handles(4),'\tau_{ORN}-\tau_{PID} (ms)','interpreter','tex')

% show the timing data for the LN model
stim_half_time = NaN(1e4,1); % time it takes to go to half max
resp_half_time = NaN(1e4,1); % time it takes to go to half max
foreground_stim = NaN(1e4,1);
c = 1;
a = 1050;
z = 1500; % nominal stimulus start and stop
for i = 1:length(PulseData)
	for j = 1:width(PulseData(i).stim)
		this_stim = PulseData(i).stim(:,j);
		foreground_stim(c) = mean(this_stim(a:z));
		this_stim = this_stim/max(this_stim(a:z));
		this_resp = PulseData(i).pred(:,j);
		this_resp = this_resp -  mean(this_resp(1:a));
		this_resp = this_resp/max(this_resp(a:z));
		stim_half_time(c) = max([find(this_stim(a:z)>.5,1,'first') NaN]);
		resp_half_time(c) = max([find(this_resp(a:z)>.5,1,'first') NaN]);
		c = c+1;
	end
end
stim_half_time(c:end) = [];
resp_half_time(c:end) = [];
foreground_stim(c:end) = [];

resp_time = resp_half_time - stim_half_time;

% clean up exactly like we did the data
foreground_stim(rm_this) = [];
resp_time(rm_this) = [];

l=plot(axes_handles(4),foreground_stim,resp_time,'r+');

% show the spearman rho with the data
s=rsquare(resp_time,resp_time_data);
set(axes_handles(4),'XScale','log','XLim',[1e-2 20],'YLim',[0 100])
legend(l,strcat('r^2=',oval(s)))

% ##       ##    ##    ##      ## ######## ########  ######## ########  
% ##       ###   ##    ##  ##  ## ##       ##     ## ##       ##     ## 
% ##       ####  ##    ##  ##  ## ##       ##     ## ##       ##     ## 
% ##       ## ## ##    ##  ##  ## ######   ########  ######   ########  
% ##       ##  ####    ##  ##  ## ##       ##     ## ##       ##   ##   
% ##       ##   ###    ##  ##  ## ##       ##     ## ##       ##    ##  
% ######## ##    ##     ###  ###  ######## ########  ######## ##     ## 

load('MSG_per_neuron.mat','MSG_data')
load('LN_Fit_to_MSG.mat','p')
% fit a LN model for each neuron
% clear p
% for i = 1:13
% 	clear d
% 	c = 1;
% 	for j = 1:8
% 		if width(MSG_data(j,i).stim)>1
% 			d(c).stimulus = mean2(MSG_data(j,i).stim);
% 			d(c).response = mean2(MSG_data(j,i).resp);
% 			c = c+1;
% 		end
% 	end
% 	if exist('d','var')
% 		p(i) = FitModel2Data(@pLNModel,d,'p0',p(i),'nsteps',2);
% 	end
% end

% compute all the model predictions and the fit qualities 
qx = NaN(8,13);
qy = NaN(8,13);
for i = 1:8 % there are 8 paradimgs in the MSG data
	for j = 1:13 % there are 13 neurons
		if ~isempty(p(j).A) && width(MSG_data(i,j).stim) > 1
			for k = 1:width(MSG_data(i,j).stim)
				MSG_data(i,j).fp_LN(:,k) = pLNModel(MSG_data(i,j).stim(:,k),p(j));
			end
			% compute the geffen-Meister error metric
			[qx(i,j) qy(i,j)] = GeffenMeister(MSG_data(i,j).resp,MSG_data(i,j).fp_LN);
		end
	end
end

% plot on the error plot, and colour code by stimulus
c = parula(9);
for i = 1:8
	this_qx = qx(i,~isnan(qx(i,:)));
	this_qy = qy(i,~isnan(qy(i,:)));
	plot(axes_handles(1),this_qx,this_qy,'x','Color',c(i,:))
end


% get the gain scaling factor -- this is to correct for the fact that we normalised the filter heights everywhere
x = mean2(horzcat(MSG_data(1,:).stim));
y = mean2(horzcat(MSG_data(1,:).resp));
gsf = std(y)/std(x);


% show the gain plot again, and also the gain for the da model 
gain = NaN(8,13);
gain_DA = NaN(8,13);
mean_stim = NaN(8,13);
for i = 1:8 % iterate over all paradigms 
	for j = 1:13
		if width(MSG_data(i,j).stim) > 1
			y = MSG_data(i,j).resp; % average over all neurons 
			x = MSG_data(i,j).fp;
			x_da = MSG_data(i,j).fp_LN;
			if ~isvector(x)
				x = mean2(x);
			end
			if ~isvector(y)
				y = mean2(y);
			end 
			if ~isvector(x_da)
				x_da = mean2(x_da);
			end 


			% trim NaNs again
			rm_this = isnan(x) | isnan(y) | isnan(x_da);
			x(rm_this) = [];
			x_da(rm_this) = [];
			y(rm_this) = [];

			temp=fit(x(:),y(:),'poly1');
			gain(i,j) = temp.p1;
			mean_stim(i,j) = mean(mean([MSG_data(i,j).stim]));

			temp=fit(x(:),x_da(:),'poly1');
			gain_DA(i,j) = temp.p1;

			% get the units of gain right
			gain(i,j) = gain(i,j)*gsf;
			gain_DA(i,j) = gain_DA(i,j)*gsf;

		end
	end	
end
plot(axes_handles(2),mean_stim(:),gain(:),'kx')
l=plot(axes_handles(2),mean_stim(:),gain_DA(:),'rx');
xlabel(axes_handles(2),'Mean Stimulus (V)')
ylabel(axes_handles(2),'Neuron Gain (Hz/V)')
set(axes_handles(2),'XScale','log','YScale','log')

% show the spearman rho 
a = gain(:);
b = gain_DA(:);
rm_this = isnan(a) | isnan(b);
a(rm_this) = [];
b(rm_this) = [];
s = rsquare(a,b);
legend(l,strcat('r^2=',oval(s)));


% ##       ##    ##     ######  ########  ######## ######## ########  ##     ## ########  
% ##       ###   ##    ##    ## ##     ## ##       ##       ##     ## ##     ## ##     ## 
% ##       ####  ##    ##       ##     ## ##       ##       ##     ## ##     ## ##     ## 
% ##       ## ## ##     ######  ########  ######   ######   ##     ## ##     ## ########  
% ##       ##  ####          ## ##        ##       ##       ##     ## ##     ## ##        
% ##       ##   ###    ##    ## ##        ##       ##       ##     ## ##     ## ##        
% ######## ##    ##     ######  ##        ######## ######## ########   #######  ##        
   


peak_loc_data = NaN(8,13);
peak_loc_DA = NaN(8,13);
for i = 1:8 % iterate over all paradigms 
	for j = 1:13
		if width(MSG_data(i,j).stim) > 1
			this_resp = MSG_data(i,j).resp;
			this_stim = MSG_data(i,j).stim;
			this_pred = MSG_data(i,j).fp_LN;
			if ~isvector(this_resp)
				this_resp = mean2(this_resp);
			end
			if ~isvector(this_stim)
				this_stim = mean2(this_stim);
			end
			if ~isvector(this_pred)
				this_pred = mean2(this_pred);
			end

			a = this_resp - mean(this_resp);
			b = this_stim - mean(this_stim);
			a = a/std(a);
			b = b/std(b);
			x = xcorr(a,b); % positive peak means a lags b
			t = 1e-3*(1:length(x));
			t = t-mean(t);
			x = x/max(x);
			[~,loc] = max(x);
			peak_loc_data(i,j) = t(loc);

			% throw out the first 5 seconds
			this_pred(1:5e3) = [];
			b(1:5e3) = [];
			a = this_pred - mean(this_pred);
			a = a/std(a);
			x = xcorr(a,b); % positive peak means a lags b
			x = x/max(x);
			t = 1e-3*(1:length(x));
			t = t-mean(t);
			[~,loc] = max(x);
			peak_loc_DA(i,j) = t(loc);
			mean_stim(i,j) = mean(mean([MSG_data(i,j).stim]));

		end
	end
end 

plot(axes_handles(3),mean_stim(:),1e3*peak_loc_data(:),'kx')
l=plot(axes_handles(3),mean_stim(:),1e3*peak_loc_DA(:),'rx');
set(axes_handles(3),'XLim',[-.1 4],'YLim',[5 120])
xlabel(axes_handles(3),'Mean Stimulus (V)')
ylabel(axes_handles(3),'Peak xcorr. (ms)')

a = peak_loc_data(:);
b = peak_loc_DA(:);
rm_this = isnan(a) | isnan(b) | a < 0 | b < 0;
a(rm_this) = [];
b(rm_this) = [];
s = rsquare(a,b);
legend(l,strcat('r^2=',oval(s)));

return



%       ########     ###               ##     ##  #######  ########  ######## ##       
%       ##     ##   ## ##              ###   ### ##     ## ##     ## ##       ##       
%       ##     ##  ##   ##             #### #### ##     ## ##     ## ##       ##       
%       ##     ## ##     ##            ## ### ## ##     ## ##     ## ######   ##       
%       ##     ## #########            ##     ## ##     ## ##     ## ##       ##       
%       ##     ## ##     ##            ##     ## ##     ## ##     ## ##       ##       
%       ########  ##     ##            ##     ##  #######  ########  ######## ######## 


% DA Model
load('LVF_data.mat')
% generate predictions
for i = 1:length(data)
	data(i).fp = DAModelv2(mean2(data(i).stim),p(i));
	% censor the first 10 seconds
	resp = data(i).resp(1e4:end,:);
	fp = data(i).fp(1e4:end);
	% show fit quality of DA model for Large Variance Flicker
	[qx, qy] = GeffenMeister(resp,fp);
	plot(axes_handles(6),qx,qy,'ko')
end


% show gain analysis -- da model
ph = []; ph(3) = axes_handles(10);
history_lengths = 0.4890;
time = 1e-3*(1:length(data(1).fp));
[p,~,~,~,~,history_lengths]=GainAnalysisWrapper2('response',mean2([data.resp]),'prediction',mean2([data.fp]),'stimulus',mean2([data.stim]),'time',time,'ph',ph,'history_lengths',history_lengths,'example_history_length',history_lengths(1),'engine',@GainAnalysis5);


% show the p-value
axes(axes_handles(10))
text(10,60,strkat('p = ',oval(p(1))))



%     ########     ###             ########  ##     ## ##        ######  ########  ######  
%     ##     ##   ## ##            ##     ## ##     ## ##       ##    ## ##       ##    ## 
%     ##     ##  ##   ##           ##     ## ##     ## ##       ##       ##       ##       
%     ##     ## ##     ##          ########  ##     ## ##        ######  ######    ######  
%     ##     ## #########          ##        ##     ## ##             ## ##             ## 
%     ##     ## ##     ##          ##        ##     ## ##       ##    ## ##       ##    ## 
%     ########  ##     ##          ##         #######  ########  ######  ########  ######  


% pulses

load('PulseData.mat')
% make all the model predictions
qx = NaN(length(PulseData),1);
qy = NaN(length(PulseData),1);
for i = 1:length(PulseData)
	for j = 1:width(PulseData(i).stim)
		PulseData(i).fp(:,j) = DAModelv2(PulseData(i).stim(:,j),p);
	end
	[qx(i), qy(i)] = GeffenMeister(PulseData(i).resp,PulseData(i).fp);
end

% show this on the plot
plot(axes_handles(6),qx,qy,'k+')


% show the timing data
stim_half_time = NaN(1e4,1); % time it takes to go to half max
resp_half_time = NaN(1e4,1); % time it takes to go to half max
foreground_stim = NaN(1e4,1);
c = 1;
a = 1050;
z = 1500; % nominal stimulus start and stop
for i = 1:length(PulseData)
	for j = 1:width(PulseData(i).stim)
		this_stim = PulseData(i).stim(:,j);
		foreground_stim(c) = mean(this_stim(a:z));
		this_stim = this_stim/max(this_stim(a:z));
		this_resp = PulseData(i).resp(:,j);
		this_resp = this_resp -  mean(this_resp(1:a));
		this_resp = this_resp/max(this_resp(a:z));
		stim_half_time(c) = max([find(this_stim(a:z)>.5,1,'first') NaN]);
		resp_half_time(c) = max([find(this_resp(a:z)>.5,1,'first') NaN]);
		c = c+1;
	end
end
stim_half_time(c:end) = [];
resp_half_time(c:end) = [];
resp_time = resp_half_time - stim_half_time;
foreground_stim(c:end) = [];

% clean up
rm_this = foreground_stim < 1e-2 | resp_time < 0 | resp_time > 100;
foreground_stim(rm_this) = [];
resp_time(rm_this) = [];
resp_time_data = resp_time;

plot(axes_handles(9),foreground_stim,resp_time,'k+')
set(axes_handles(9),'XScale','log','XLim',[1e-2 20],'YLim',[0 100])
xlabel(axes_handles(9),'Mean Stimulus (V)')
ylabel(axes_handles(9),'\tau_{ORN}-\tau_{PID} (ms)','interpreter','tex')

% show the timing data for the DA model
stim_half_time = NaN(1e4,1); % time it takes to go to half max
resp_half_time = NaN(1e4,1); % time it takes to go to half max
foreground_stim = NaN(1e4,1);
c = 1;
a = 1050;
z = 1500; % nominal stimulus start and stop
for i = 1:length(PulseData)
	for j = 1:width(PulseData(i).stim)
		this_stim = PulseData(i).stim(:,j);
		foreground_stim(c) = mean(this_stim(a:z));
		this_stim = this_stim/max(this_stim(a:z));
		this_resp = PulseData(i).fp(:,j);
		this_resp = this_resp -  mean(this_resp(1:a));
		this_resp = this_resp/max(this_resp(a:z));
		stim_half_time(c) = max([find(this_stim(a:z)>.5,1,'first') NaN]);
		resp_half_time(c) = max([find(this_resp(a:z)>.5,1,'first') NaN]);
		c = c+1;
	end
end
stim_half_time(c:end) = [];
resp_half_time(c:end) = [];
foreground_stim(c:end) = [];

resp_time = resp_half_time - stim_half_time;

% clean up exactly like we did the data
foreground_stim(rm_this) = [];
resp_time(rm_this) = [];

l=plot(axes_handles(9),foreground_stim,resp_time,'r+');

s=rsquare(resp_time,resp_time_data);
set(axes_handles(4),'XScale','log','XLim',[1e-2 20],'YLim',[0 100])
legend(l,strcat('r^2=',oval(s)))


%      ########     ###              ##      ## ######## ########  ######## ########  
%      ##     ##   ## ##             ##  ##  ## ##       ##     ## ##       ##     ## 
%      ##     ##  ##   ##            ##  ##  ## ##       ##     ## ##       ##     ## 
%      ##     ## ##     ##           ##  ##  ## ######   ########  ######   ########  
%      ##     ## #########           ##  ##  ## ##       ##     ## ##       ##   ##   
%      ##     ## ##     ##           ##  ##  ## ##       ##     ## ##       ##    ##  
%      ########  ##     ##            ###  ###  ######## ########  ######## ##     ## 

% show weber scaling 
load('MSG_per_neuron.mat','MSG_data')
load('DA_Fit_to_MSG.mat')

% compute all the model predictions and the fit qualities 
qx = NaN(8,13);
qy = NaN(8,13);
for i = 1:8 % there are 8 paradimgs in the MSG data
	for j = 1:13 % there are 13 neurons
		if ~isempty(p(j).A) && width(MSG_data(i,j).stim) > 1
			for k = 1:width(MSG_data(i,j).stim)
				MSG_data(i,j).fp_DA(:,k) = DAModelv2(MSG_data(i,j).stim(:,k),p(j));
			end
			% compute the geffen-Meister error metric
			[qx(i,j) qy(i,j)] = GeffenMeister(MSG_data(i,j).resp,MSG_data(i,j).fp_DA);
		end
	end
end

% plot on the error plot, and colour code by stimulus
c = parula(9);
for i = 1:8
	this_qx = qx(i,~isnan(qx(i,:)));
	this_qy = qy(i,~isnan(qy(i,:)));
	plot(axes_handles(6),this_qx,this_qy,'x','Color',c(i,:))
end


% get the gain scaling factor -- this is to correct for the fact that we normalised the filter heights everywhere
x = mean2(horzcat(MSG_data(1,:).stim));
y = mean2(horzcat(MSG_data(1,:).resp));
gsf = std(y)/std(x);


% show the gain plot again, and also the gain for the da model 
gain = NaN(8,13);
gain_DA = NaN(8,13);
mean_stim = NaN(8,13);
for i = 1:8 % iterate over all paradigms 
	for j = 1:13
		if width(MSG_data(i,j).stim) > 1
			y = MSG_data(i,j).resp; % average over all neurons 
			x = MSG_data(i,j).fp;
			x_da = MSG_data(i,j).fp_DA;
			if ~isvector(x)
				x = mean2(x);
			end
			if ~isvector(y)
				y = mean2(y);
			end 
			if ~isvector(x_da)
				x_da = mean2(x_da);
			end 


			% trim NaNs again
			rm_this = isnan(x) | isnan(y) | isnan(x_da);
			x(rm_this) = [];
			x_da(rm_this) = [];
			y(rm_this) = [];

			temp=fit(x(:),y(:),'poly1');
			gain(i,j) = temp.p1;
			mean_stim(i,j) = mean(mean([MSG_data(i,j).stim]));

			temp=fit(x(:),x_da(:),'poly1');
			gain_DA(i,j) = temp.p1;

			% get the units of gain right
			gain(i,j) = gain(i,j)*gsf;
			gain_DA(i,j) = gain_DA(i,j)*gsf;

		end
	end	
end
plot(axes_handles(7),mean_stim(:),gain(:),'kx')
l=plot(axes_handles(7),mean_stim(:),gain_DA(:),'rx');
xlabel(axes_handles(7),'Mean Stimulus (V)')
ylabel(axes_handles(7),'Neuron Gain (Hz/V)')
set(axes_handles(7),'XScale','log','YScale','log')

% show the spearman rho 
a = gain(:);
b = gain_DA(:);
rm_this = isnan(a) | isnan(b);
a(rm_this) = [];
b(rm_this) = [];
s = rsquare(a,b);
legend(l,strcat('r^2=',oval(s)));




%    ########     ###        ######  ########  ######## ######## ########  ##     ## ########  
%    ##     ##   ## ##      ##    ## ##     ## ##       ##       ##     ## ##     ## ##     ## 
%    ##     ##  ##   ##     ##       ##     ## ##       ##       ##     ## ##     ## ##     ## 
%    ##     ## ##     ##     ######  ########  ######   ######   ##     ## ##     ## ########  
%    ##     ## #########          ## ##        ##       ##       ##     ## ##     ## ##        
%    ##     ## ##     ##    ##    ## ##        ##       ##       ##     ## ##     ## ##        
%    ########  ##     ##     ######  ##        ######## ######## ########   #######  ##        


% xcorr for mean shifted gaussians
peak_loc_data = NaN(8,13);
peak_loc_DA = NaN(8,13);
for i = 1:8 % iterate over all paradigms 
	for j = 1:13
		if width(MSG_data(i,j).stim) > 1
			this_resp = MSG_data(i,j).resp;
			this_stim = MSG_data(i,j).stim;
			this_pred = MSG_data(i,j).fp_DA;
			if ~isvector(this_resp)
				this_resp = mean2(this_resp);
			end
			if ~isvector(this_stim)
				this_stim = mean2(this_stim);
			end
			if ~isvector(this_pred)
				this_pred = mean2(this_pred);
			end

			a = this_resp - mean(this_resp);
			b = this_stim - mean(this_stim);
			a = a/std(a);
			b = b/std(b);
			x = xcorr(a,b); % positive peak means a lags b
			t = 1e-3*(1:length(x));
			t = t-mean(t);
			x = x/max(x);
			[~,loc] = max(x);
			peak_loc_data(i,j) = t(loc);

			% throw out the first 5 seconds
			this_pred(1:5e3) = [];
			b(1:5e3) = [];
			a = this_pred - mean(this_pred);
			a = a/std(a);
			x = xcorr(a,b); % positive peak means a lags b
			x = x/max(x);
			t = 1e-3*(1:length(x));
			t = t-mean(t);
			[~,loc] = max(x);
			peak_loc_DA(i,j) = t(loc);
			mean_stim(i,j) = mean(mean([MSG_data(i,j).stim]));

		end
	end
end 

plot(axes_handles(8),mean_stim(:),1e3*peak_loc_data(:),'kx')
l=plot(axes_handles(8),mean_stim(:),1e3*peak_loc_DA(:),'rx');
set(axes_handles(8),'XLim',[-.1 4],'YLim',[5 120])
xlabel(axes_handles(8),'Mean Stimulus (V)')
ylabel(axes_handles(8),'Peak xcorr. (ms)')

a = peak_loc_data(:);
b = peak_loc_DA(:);
rm_this = isnan(a) | isnan(b) | a < 0 | b < 0;
a(rm_this) = [];
b(rm_this) = [];
s = rsquare(a,b);
legend(l,strcat('r^2=',oval(s)));


prettyFig('plw=1.5;','lw=1.5;','fs=14;')

if being_published
	snapnow
	delete(gcf)
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

t = toc;
%% 
% This document was built in: 
disp(strcat(oval(t,3),' seconds.'))



