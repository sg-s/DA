% Explaining_Weber_Gain.m
% 
% created by Srinivas Gorur-Shandilya at 2:16 , 22 January 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


pHeader;

%% Explaining ORN gain changes to Mean Shifted Gaussians
% In this document, we look at how the gain of a ORN changes if we change the mean of the stimulus, and try to understand why it changes in the way it does. 

% get the data and prepare it


[PID, ~, fA, paradigm, orn, fly, AllControlParadigms, paradigm_hashes] = consolidateData('/local-data/DA-paper/LFP-MSG/september',1);

% remove baseline from all PIDs
for i = 1:width(PID)
	PID(:,i) = PID(:,i) - mean(PID(1:5e3,i));
end

% sort the paradigms sensibly
sort_value = [];
for i = 1:length(AllControlParadigms)
	sort_value(i) = (mean(AllControlParadigms(i).Outputs(1,:)));
end
[~,idx] = sort(sort_value);


AllControlParadigms = AllControlParadigms(idx);
paradigm_new = paradigm*NaN;
for i = 1:length(idx)
	paradigm_new(paradigm == idx(i)) = i;
end
paradigm = paradigm_new;

% throw our bad traces
bad_trials = (sum(fA) == 0 | isnan(sum(fA)));
PID(:,bad_trials) = [];
fA(:,bad_trials) = [];
paradigm(bad_trials) = [];
orn(bad_trials) = [];

% some core variables
c = parula(max(paradigm)+1); % colour scheme
dt = 1e-3;
time = dt*(1:length(PID));

% extract filters and find gain
a = 35e3; z = 55e3;
[K,fA_pred,fA_gain] = extractFilters(PID,fA,'use_cache',true,'a',a,'z',z);

% use a single filter to re-project the stimulus
K0 = mean(K(:,paradigm==1),2);
filtertime = dt*(1:length(K0)) - .1;
proj_stim = NaN*fA_pred; 
for i = 1:width(PID)
	proj_stim(:,i) = convolve(time,PID(:,i),K0,filtertime);
end


% show gain changes for all paradigms -- average over neurons 
ss = 100;
all_x = 0:0.1:2;
for i = 1:max(paradigm) % iterate over all paradigms 
	y = nanmean(fA(a:z,paradigm == i),2);
	x = nanmean(fA_pred(a:z,paradigm == i),2);
	s = nanmean(proj_stim(a:z,paradigm == i),2);
	x = x - nanmean(x);
	x = x + nanmean(nanmean(s));
	[~,orn_io_data(i)] = plotPieceWiseLinear(x,y,'nbins',50,'Color',c(i,:),'make_plot',false);
end

%% Ansatz: Matching Stimulus Statistics
% In this section, we assume that the ORNs gain changes in an adaptive manner, to better encode the stimulus. 


figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[800 800]); hold on
subplot(3,3,1), hold on
x = 0:0.05:2;
y = zeros(length(x),max(paradigm));
for i = 1:max(paradigm)
	stim = proj_stim(a:z,paradigm == i);
	[hy,hx]  = hist(stim(:),50);
	temp = histcounts(stim(:),x);
	temp = temp/sum(temp);
	temp = cumsum(temp);
	y(2:end,i) = temp;
	hy = hy/sum(hy);
	plot(hx,hy,'Color',c(i,:));
end
ylabel('p(Proj. stimulus)')
set(gca,'XLim',[0 1.5])
title('Stimulus-defined I/O curves')

% normalise ORN data
m = max(max([orn_io_data.y]));
for i = 1:max(paradigm)
	orn_io_data(i).y = orn_io_data(i).y/m;
end

% now plot the integrals -- these are the theoretical predictions
% also plot the data
subplot(3,3,4), hold on
for i = 1:max(paradigm)
	plot(x,y(:,i),'Color',c(i,:));
	plot(orn_io_data(i).x,orn_io_data(i).y,'+','Color',c(i,:))
end
xlabel('Proj. Stimulus (V)')
ylabel('Response (norm)')
set(gca,'XLim',[0 1.5])

% now show the correlation between predicted gain and actual gain
subplot(3,3,7), hold on
orn_gain = NaN(max(paradigm),1);
predicted_gain = NaN*orn_gain;
for i = 1:max(paradigm)
	temp = interp1(orn_io_data(i).x,orn_io_data(i).y,x);
	temp = fit(x(~isnan(temp))',temp(~isnan(temp))','poly1');
	orn_gain(i) = temp.p1;
	temp = interp1(orn_io_data(i).x,orn_io_data(i).y,x);
	temp = fit(x(~isnan(temp))',y(~isnan(temp),i),'poly1');
	predicted_gain(i) = temp.p1;
end
plot(predicted_gain,orn_gain,'k+')
plot([1e-3 10],[1e-3 10],'k--')
set(gca,'XScale','log','YScale','log','YLim',[.1 10])
xlabel('Predicted Gain')
ylabel('ORN Gain')

% now let Kd and n vary, and try to match the theoretical stimulus-defined distributions.
clear p
p = cache('unconstrained_MSG_fit2');
if isempty(p)
	clear d
	d.stimulus = x;
	d.response = y(:,1)';
	p = fitModel2Data(@hill,d,'make_plot',false,'nsteps',1000,'display_type','iter');
	for i = 2:max(paradigm)
		clear d
		d.stimulus = x;
		d.response = y(:,i)';
		p(i) = fitModel2Data(@hill,d,'p0',p(1),'make_plot',false,'nsteps',300,'display_type','iter');
	end
	cache('unconstrained_MSG_fit2',p);
end

subplot(3,3,2), hold on

for i = 1:max(paradigm)
	plot(x,y(:,i),'Color',c(i,:));
	pred = hill(x,p(i));
	plot(x,pred,'LineStyle','--','Color',c(i,:))
end
title('K_D and n can vary')
set(gca,'XLim',[0 2])
xlabel('Stimulus (V)')
ylabel('Response (norm)')


subplot(3,3,5), hold on
for i = 1:max(paradigm)
	pred = hill(x,p(i));
	plot(x,pred,'LineStyle','--','Color',c(i,:))
	plot(orn_io_data(i).x,orn_io_data(i).y,'+','Color',c(i,:))
end
set(gca,'XLim',[0 2])
xlabel('Stimulus (V)')
ylabel('Response (norm)')

% now show the correlation between predicted gain and actual gain
subplot(3,3,8), hold on
orn_gain = NaN(max(paradigm),1);
predicted_gain = NaN*orn_gain;
for i = 1:max(paradigm)
	temp = interp1(orn_io_data(i).x,orn_io_data(i).y,x);
	temp = fit(x(~isnan(temp))',temp(~isnan(temp))','poly1');
	orn_gain(i) = temp.p1;
	temp = interp1(orn_io_data(i).x,orn_io_data(i).y,x);
	pred = hill(x,p(i));
	temp = fit(x(~isnan(temp))',pred(~isnan(temp))','poly1');
	predicted_gain(i) = temp.p1;
end
plot(predicted_gain,orn_gain,'k+')
plot([1e-3 10],[1e-3 10],'k--')
set(gca,'XScale','log','YScale','log','YLim',[.1 10])
xlabel('Predicted Gain')


% now plot the stimulus + best fit assuming Kd constraint
clear p
p = cache('Kd_constrained_MSG_fit2');
if isempty(p)
	clear d
	d.stimulus = x;
	d.response = y(:,1)';
	p = fitModel2Data(@hill,d,'make_plot',false,'nsteps',1000,'display_type','iter');
	ub = p;
	ub.k = Inf;
	lb = p;
	lb.k = 0;
	for i = 2:max(paradigm)
		clear d
		d.stimulus = x;
		d.response = y(:,i)';
		p(i) = fitModel2Data(@hill,d,'p0',p(1),'ub',ub,'lb',lb,'make_plot',false,'nsteps',300,'display_type','iter');
	end
	cache('Kd_constrained_MSG_fit2',p);
end

subplot(3,3,3), hold on
for i = 1:max(paradigm)
	plot(x,y(:,i),'Color',c(i,:));
	pred = hill(x,p(i));
	plot(x,pred,'LineStyle','--','Color',c(i,:))
end
title('Only K_D can vary')
set(gca,'XLim',[0 2])
xlabel('Stimulus (V)')
ylabel('Response (norm)')

subplot(3,3,6), hold on
for i = 1:max(paradigm)
	pred = hill(x,p(i));
	plot(x,pred,'LineStyle','--','Color',c(i,:))
	plot(orn_io_data(i).x,orn_io_data(i).y,'+','Color',c(i,:))
end
set(gca,'XLim',[0 2])
xlabel('Stimulus (V)')
ylabel('Response (norm)')

% now show the correlation between predicted gain and actual gain
subplot(3,3,9), hold on
orn_gain = NaN(max(paradigm),1);
predicted_gain = NaN*orn_gain;
for i = 1:max(paradigm)
	temp = interp1(orn_io_data(i).x,orn_io_data(i).y,x);
	temp = fit(x(~isnan(temp))',temp(~isnan(temp))','poly1');
	orn_gain(i) = temp.p1;
	temp = interp1(orn_io_data(i).x,orn_io_data(i).y,x);
	pred = hill(x,p(i));
	temp = fit(x(~isnan(temp))',pred(~isnan(temp))','poly1');
	predicted_gain(i) = temp.p1;
end
plot(predicted_gain,orn_gain,'k+')
plot([1e-3 10],[1e-3 10],'k--')
set(gca,'XScale','log','YScale','log','YLim',[.1 10])
xlabel('Predicted Gain')

prettyFig('plw=1.3;','lw=1.5;','fs=12;')

if being_published
	snapnow
	delete(gcf)
end



%% Version Info
%
pFooter;


