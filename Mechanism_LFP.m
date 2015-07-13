% Mechanism_LFP
% attempt to understand mechanism of slow adaptation using the LFP
% 
% created by Srinivas Gorur-Shandilya at 4:04 , 10 July 2015. Contact me at http://srinivas.gs/contact/
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

%% What does the LFP tell us about the mechanism of slow adaptation?
% In this document, we attempt to figure out something about the mechanism of slow adaptation by looking at how the LFP and firing rates change in response to fluctuating odour stimuli with different means but similar variances. 


%      ######  ######## #### ##     ## ##     ## ##       ##     ##  ######  
%     ##    ##    ##     ##  ###   ### ##     ## ##       ##     ## ##    ## 
%     ##          ##     ##  #### #### ##     ## ##       ##     ## ##       
%      ######     ##     ##  ## ### ## ##     ## ##       ##     ##  ######  
%           ##    ##     ##  ##     ## ##     ## ##       ##     ##       ## 
%     ##    ##    ##     ##  ##     ## ##     ## ##       ##     ## ##    ## 
%      ######     ##    #### ##     ##  #######  ########  #######   ######  


%% Stimulus Characterisation
% First, we show that we are able to deliver Gaussian-distributed odour stimuli, and that we are able to vary the means of the distributions of these stimuli. 


[PID, LFP, fA, paradigm, orn, AllControlParadigms, paradigm_hashes] = consolidateData('/local-data/DA-paper/LFP-MSG/fig1-flicker/',1);


% sort the paradigms sensibly
sort_value = [];
for i = 1:length(AllControlParadigms)
	sort_value(i) = mean(mean(AllControlParadigms(i).Outputs));
end

[~,idx] = sort(sort_value);
AllControlParadigms = AllControlParadigms(idx);
paradigm_new = paradigm*NaN;
for i = 1:length(idx)
	paradigm_new(paradigm == idx(i)) = i;
end
paradigm = paradigm_new;


% throw out trials where we didn't record the LFP, for whatever reason
not_LFP = 0*orn;
for i = 1:width(LFP)
	not_LFP(i) = abs(mean2(LFP(:,i)));
end
LFP(:,not_LFP< 0.5) = NaN;


% remove linear trend for all PID
for i = 1:width(PID)
	y = PID(30e3:55e3,i); 
	x = 1:length(y); x = x(:);
	ff = fit(x,y,'poly1');
	PID(30e3:55e3,i) = y - ff(x) + mean(y);
end

%%
% The following figure shows the distribution of the inputs, for the terminal 20seconds of the flickering odour stimulus. 



figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
c = parula(1+length(unique(paradigm)));
for i = 1:length(unique(paradigm))
	hist_this = PID(30e3:55e3,paradigm==i);
	xx =  linspace(min(min(hist_this)),max(max(hist_this)),50);
	y = NaN(sum(paradigm==i),50);
	for j = 1:sum(paradigm==i)
		y(j,:) = hist(hist_this(:,j),xx);
		y(j,:) = y(j,:)/sum(y(j,:));
	end
	% get everything on the same x axis
	if width(y) > 1
		errorShade(xx,mean(y),std(y)/sqrt(width(y)),'Color',c(i,:),'LineWidth',5);
	else
		plot(xx,y,'Color',c(i,:));
	end
end
xlabel('Stimulus (V)')
ylabel('p(stimulus)')
PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%%
% How reproducible is the stimulus? In the following figure, we show the same stimulus from the entire dataset, over all the trails, and all the flies we have. The shading shows the standard error of the mean.

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot_this = PID(40e3:end,paradigm==1);

time = 1e-3*(1:length(plot_this)) + 40;
errorShade(time,mean2(plot_this),std(plot_this')/sqrt(width(plot_this)));
xlabel('Time (s)')
ylabel('PID (V)')
title(strcat('Stimulus reproducibility: n = ',oval(width(plot_this))))

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end


% 								##       ######## ########  
% 								##       ##       ##     ## 
% 								##       ##       ##     ## 
% 								##       ######   ########  
% 								##       ##       ##        
% 								##       ##       ##        
% 								######## ##       ##        


%% Local Field Potential
% We now look at the LFP. Here is an example neuron, showing how the LFP changes with the different paradigms. In the following figure, we plot the raw LFP traces, downsampled to 1kHz from the actual 10kHz trace, band passed filtered to remove slow fluctuations and spikes. It is immediately clear that the magnitude of fluctuations of the LFP goes down with increasing mean of the stimulus. 

example_orn = 1;
figure('outerposition',[0 0 1200 500],'PaperUnits','points','PaperSize',[1200 500]); hold on
these_paradigms = unique(paradigm(:,orn == example_orn));
c = parula(1+length(these_paradigms));
time = 1e-3*(1:length(LFP(30e3:40e3,1))) + 30;
l = [];
L = {};
for i = 1:length(these_paradigms)
	plot_this = find(orn == example_orn & paradigm == these_paradigms(i));
	for j = 1:length(plot_this)
		this_LFP = LFP(:,plot_this(j));
		this_LFP = this_LFP(30e3:40e3);
		this_LFP = filter_trace(this_LFP,1000,10);
		l(i) = plot(time,this_LFP,'Color',c(i,:));
		L{i} = AllControlParadigms(i).Name;
	end
end

xlabel('Time (s)')
ylabel('\DeltaLFP (100x V)')
legend(l,L,'Location','southeastoutside')
PrettyFig;

if being_published
	snapnow
	delete(gcf)
end


%%
% We can clearly see that the fluctuations in the LFP decrease with increasing stimulus mean, suggesting decreasing gain. 

%            ######## #### ##       ######## ######## ########   ######  
%            ##        ##  ##          ##    ##       ##     ## ##    ## 
%            ##        ##  ##          ##    ##       ##     ## ##       
%            ######    ##  ##          ##    ######   ########   ######  
%            ##        ##  ##          ##    ##       ##   ##         ## 
%            ##        ##  ##          ##    ##       ##    ##  ##    ## 
%            ##       #### ########    ##    ######## ##     ##  ######  


%%
% To see this more clearly, we back out all possible filters from the stimulus, LFP and the ORN response. In the following figure, filters for the first ORN are plotted: 

if ~exist('K','var')
	K = NaN(1e3,length(orn));
	for i = 1:length(orn)
		resp = LFP(30e3:55e3,i);
		rm_this = isnan(resp);
		resp(rm_this) = [];
		if length(resp) 
			resp = filter_trace(resp,1e3,10);
			stim = PID(30e3:55e3,i);
			stim(rm_this) = [];
			stim(1:500) = [];
			resp(end-499:end) = [];
			K(:,i) = FitFilter2Data(stim,resp,[],'reg=1;','filter_length=999;');
		end
	end
end

if ~exist('K2','var')
	K2 = NaN(1e3,length(orn));
	for i = 1:length(orn)
		resp = fA(30e3:55e3,i);
		stim = LFP(30e3:55e3,i);
		rm_this = isnan(resp) | isnan(stim);
		resp(rm_this) = [];
		stim(rm_this) = [];
		if length(resp) 
			stim = filter_trace(stim,1e3,10);
			stim(1:500) = [];
			resp(end-499:end) = [];
			K2(:,i) = FitFilter2Data(stim,resp,[],'reg=1;','filter_length=999;');
		end
	end
end

if ~exist('K3','var')
	K3 = NaN(1e3,length(orn));
	for i = 1:length(orn)
		resp = fA(30e3:55e3,i);
		stim = PID(30e3:55e3,i);
		rm_this = isnan(resp) | isnan(stim);
		resp(rm_this) = [];
		stim(rm_this) = [];
		if length(resp) & sum(resp) > 0
			stim(1:500) = [];
			resp(end-499:end) = [];
			K3(:,i) = FitFilter2Data(stim,resp,[],'reg=1;','filter_length=999;');
		end
	end
end

c= parula(max(paradigm)+1);
l = [];
figure('outerposition',[0 0 1400 450],'PaperUnits','points','PaperSize',[1400 450]); hold on
subplot(1,3,1), hold on
title(['Neuron ' mat2str(1)])
for i = 1:max(paradigm)
	time = 1e-3*(1:501)-.2;
	plot_this = find(paradigm == i & orn == 1);
	plot_this = setdiff(plot_this,find(isnan(sum(K))));
	if length(plot_this) > 1
		l(i) = errorShade(time,mean2(K(300:800,plot_this)),std(K(300:800,plot_this)')/length(plot_this),'Color',c(i,:));
	elseif length(plot_this) > 0
		l(i) = plot(time,K(300:800,plot_this),'Color',c(i,:));
	end
end
xlabel('Lag (s)')
ylabel('PID \rightarrow LFP Filter')
subplot(1,3,2), hold on
title(['Neuron ' mat2str(1)])
for i = 1:max(paradigm)
	time = 1e-3*(1:501)-.2;
	plot_this = find(paradigm == i & orn == 1);
	plot_this = setdiff(plot_this,find(isnan(sum(K2))));
	if length(plot_this) > 1
		l(i) = errorShade(time,mean2(K2(300:800,plot_this)),std(K2(300:800,plot_this)')/length(plot_this),'Color',c(i,:));
	elseif length(plot_this) > 0
		l(i) = plot(time,K2(300:800,plot_this),'Color',c(i,:));
	end
end
xlabel('Lag (s)')
ylabel('LFP \rightarrow Firing rate Filter')
subplot(1,3,3), hold on
title(['Neuron ' mat2str(1)])
for i = 1:max(paradigm)
	time = 1e-3*(1:501)-.2;
	plot_this = find(paradigm == i & orn == 1);
	plot_this = setdiff(plot_this,find(isnan(sum(K3))));
	if length(plot_this) > 1
		l(i) = errorShade(time,mean2(K3(300:800,plot_this)),std(K3(300:800,plot_this)')/length(plot_this),'Color',c(i,:));
	elseif length(plot_this) > 0
		l(i) = plot(time,K3(300:800,plot_this),'Color',c(i,:));
	end
end
xlabel('Lag (s)')
ylabel('PID \rightarrow Firing rate Filter')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end



%%
% In the following figure, we calculate the gain of the LFP using standard methods (as in Fig 1) and plot the gain vs. the mean stimulus on a trial-by-trial basis, for all the neurons in the dataset. 

mean_stim = NaN*orn;
LFP_filter_height = NaN*orn;
LFP_gain = NaN*orn;
ORN_gain = NaN*orn;

for i = 1:length(orn)
	mean_stim(i) = mean(PID(30e3:55e3,i));
	LFP_filter_height(i) = -min(K(400:600,i));
	resp = LFP(30e3:55e3,i);
	stim = PID(30e3:55e3,i);
	rm_this = isnan(resp);
	resp(rm_this) = [];
	stim(rm_this) = [];
	if length(resp)
		resp = filter_trace(resp,1e3,10);
		% make prediction
		time = 1e-3*(1:length(stim));
		this_K = K(300:800,i);
		this_K = -this_K/min(this_K);
		filtertime = 1e-3*(1:length(this_K)) - .2;
		pred = convolve(time,stim,this_K,filtertime);
		rm_this = isnan(pred);
		pred(rm_this) = []; resp(rm_this) = [];
		ff = fit(pred(:),resp(:),'poly1');
		LFP_gain(i) =  ff.p1;
	end

	resp = fA(30e3:55e3,i);
	stim = PID(30e3:55e3,i);
	rm_this = isnan(resp);
	resp(rm_this) = [];
	stim(rm_this) = [];
	if length(resp) & sum(resp)>0
		% make prediction
		time = 1e-3*(1:length(stim));
		this_K = K3(300:800,i);
		this_K = this_K/max(this_K);
		filtertime = 1e-3*(1:length(this_K)) - .2;
		pred = convolve(time,stim,this_K,filtertime);
		rm_this = isnan(pred);
		pred(rm_this) = []; resp(rm_this) = [];
		ff = fit(pred(:),resp(:),'poly1');
		ORN_gain(i) =  ff.p1;
	end

end



figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
for i = 1:max(orn)
	plot(mean_stim(orn==i),LFP_gain(orn==i),'k+')
end

clear l
cf = fit(mean_stim(:),LFP_gain(:),'power1');
l(1)=plot(gca,sort(mean_stim),cf(sort(mean_stim)),'r');
L = strcat('\beta=',oval(cf.b));
xlabel('Mean Stimulus (V)')
ylabel('LFP Gain')
set(gca,'XScale','linear','YScale','linear')
legend(l,L)

subplot(1,2,2), hold on
for i = 1:max(orn)
	plot(mean_stim(orn==i),ORN_gain(orn==i),'k+')
end
cf = fit(mean_stim(~isnan(ORN_gain))',ORN_gain(~isnan(ORN_gain))','power1');
l(1)=plot(gca,sort(mean_stim),cf(sort(mean_stim)),'r');
L = strcat('\beta=',oval(cf.b));
xlabel('Mean Stimulus (V)')
ylabel('ORN Gain')
set(gca,'XScale','linear','YScale','linear')
legend(l,L)

PrettyFig;

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
	eval(strjoin({'tag -a published',which(mfilename)}))
end
