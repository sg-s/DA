% StepOnResponses.m
% Analyses the responses of ORNs to on-steps of odor
% 
% created by Srinivas Gorur-Shandilya at 11:18 , 20 November 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

%% Analysis of Long Steps
% In this document, we analyse the response of ORNs to long on-steps of odor. The data comes from many sources. The first bit of data comes from the mean shifted gaussian experiment, and we use in initial few seconds of data from the experiment. 

calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end
% redo = 1; deliberately unset

%% Mean Shifted Gaussians
% What do the responses of the ORN to the step ON of the odor look like? The MeanShiftedGaussians dataset comes from two different experiments: neurons 1-4 are exposed to concentrations 3%, 5.5% and 6.25%. Neurons 5-8 are exposed to 3%, 3.5%, 4% and and 4.5%. The following figure shows the stimulus on the left and the response on the right, for the various heights of the step ON, grouped by these experiments. 

load('MeanShiftedGaussians.mat')


% add the experiment sub field
combined_data.experiment = NaN*combined_data.neuron;
combined_data.experiment(combined_data.neuron > 8) = 3;
combined_data.experiment(combined_data.neuron > 4 & combined_data.neuron < 9) = 2;
combined_data.experiment(combined_data.neuron < 5) = 1;

dt = 1e-3;

a = floor(4/dt);
z = floor(12/dt);
z1 = floor(11/dt);
sdata = struct(); % summary data
sdata(3).max_pid = [];
sdata(3).max_f = [];
sdata(3).end_pid = [];
sdata(3).end_f = [];


figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1), hold on
c = parula(length(paradigm_names));
for i = 1:length(paradigm_names)

	for j = 1:max(combined_data.experiment)
		plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
		plot_these = intersect(plot_these, find(combined_data.experiment==j));


		if ~isempty(plot_these)
			this_pid=mean2(combined_data.PID(plot_these,:));
			time = dt*(1:length(this_pid));
			sdata(j).end_pid = [sdata(j).end_pid mean(this_pid(z1:z))];
			time = time(a:z);
			this_pid = this_pid(a:z);
			sdata(j).max_pid = [sdata(j).max_pid max(this_pid)];
			plot(time,this_pid,'Color',c(i,:))
		end

	end
end
xlabel('Time (s)')
ylabel('Stimulus (V)')


subplot(1,3,2), hold on
for i = 1:length(paradigm_names)
	for j = 1:max(combined_data.experiment)
		plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
		plot_these = intersect(plot_these, find(combined_data.experiment==j));


		if ~isempty(plot_these)
			this_orn=mean2(combined_data.fA(:,plot_these));
			time = dt*(1:length(this_orn));
			sdata(j).end_f = [sdata(j).end_f mean(this_orn(z1:z))];
			time = time(a:z);
			this_orn = this_orn(a:z);
			sdata(j).max_f = [sdata(j).max_f max(this_orn)];
			plot(time,this_orn,'Color',c(i,:))
		end

	end

end
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')


subplot(1,3,3), hold on
for j = 1:max(combined_data.experiment)
	x = sdata(j).max_pid;
	[sx,so] = sort(x);
	plot(sx,sdata(j).max_f(so),'+-r')
	x = sdata(j).end_pid;
	[sx,so] = sort(x);
	plot(sx,sdata(j).end_f(so),'+-k')
end

xlabel('Stimulus (V)')
ylabel('Firing Rate (Hz)')
set(gca,'YLim',[0 150])

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end


%%
% What are the timescales of relaxation of firing rates in each of these cases? In the following figure, we fit a 2-parameter exponential to the firing rates from peak onwards, not to claim that the decay is exponential, but to extract a timescale for each case. In the following plots, the color of the data indicates the experiment it comes from, while the fits are shown in red. The data is grouped by neuron; so the same set of neurons appears in all data coloured the same color.  


F_adap = []; % degree of adaptation
tau_adap = []; % result of fitting 1 exp
tau_1 = []; % result of fitting 2 exp -- first term
tau_2 = []; 
experiment = [];


clear ff lh L
L = {};
figure('outerposition',[0 0 1500 800],'PaperUnits','points','PaperSize',[1500 800]); hold on
c = parula(max(combined_data.experiment));

for i = 1:length(paradigm_names)
	lh = []; L = {};
	subplot(2,4,i), hold on
	for j = 1:max(combined_data.experiment)
		plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
		plot_these = intersect(plot_these, find(combined_data.experiment==j));

		if ~isempty(plot_these)
			this_resp=mean2(combined_data.fA(:,plot_these));
			time = dt*(1:length(this_resp));
			time = time(a:z);
			this_resp = this_resp(a:z);
			[~,loc]=max(this_resp);
			t = time(loc:end);
			t = t -min(t);
			ft=fittype('a*exp(-x./b1)+ a*exp(-x./b2) +c');
			fo = fitoptions(ft);
			fo.Robust = 'on';
			tf = mean(this_resp(end-100:end));
			pf = max(this_resp);  pf = pf - tf; pf = pf/2;
			fo.Lower = [pf-1e-3 1e-2 1e-2 tf-1e-3];
			fo.Upper = [pf+1e-3 2 20 tf+1e-3];
			fo.StartPoint = [pf .1 1 tf];
			ff = fit(t(:),this_resp(loc:end),ft,fo);
			plot(time,this_resp,'Color',c(j,:))
			lh = [lh plot(time(loc:end),ff(t),'r')];
			L = [L (strcat('\tau=',oval(ff.b1,1),',',oval(ff.b2,1),'s'))];

			F_adap = [F_adap (max(this_resp) - mean(this_resp(end-1000:end)))/(max(this_resp))];
			tau_1 = [tau_1 ff.b1];
			tau_2 = [tau_2 ff.b2];

			% also fit a single term exp
			ft=fittype('a*exp(-x./b1) +c');
			fo = fitoptions(ft);
			fo.Robust = 'on';
			fo.StartPoint = [pf .1 tf];
			ff = fit(t(:),this_resp(loc:end),ft,fo);
			tau_adap = [tau_adap ff.b1];
			experiment = [experiment j];


		end
	end
	legend(lh,L)
	title(paradigm_names{i}(strfind(paradigm_names{i},'-')+1:end))
	xlabel('Time (s)')
	ylabel('Firing Rate (Hz)')
	set(gca,'YLim',[0 120])

end

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%% The Liu-Wang Model
% The Liu-Wang model predicts that the adaptation degree (ratio of peak firing - ss firing to peak firing) depends linearly on the timescale of adaptation, and that the slope of this line is the negative timescale of calcium in the ORN. Here, we fit this data to the Liu-Wang prediction:

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(tau_adap,F_adap,'k+')
ff = fit(tau_adap(:),F_adap(:),'poly1');
plot(tau_adap,ff(tau_adap),'r')
legend(strcat('\tau_{Ca}=',oval(-ff.p1*1000),'ms'))
xlabel('Adaptation Time (s)')
ylabel('Degree of adaptation')
PrettyFig;

if being_published
	snapnow
	delete(gcf)
end


%%
% These different timescales, and the dependence of timescale on the odor concentration, raises the question of whether these curves "rescale" á la Martelli et al. In the following figure, we rescale all these plots by the peak and plot them one on top of the other. In the following plots, we group the data so that all curves come from the same neuron. 

figure('outerposition',[0 0 1400 500],'PaperUnits','points','PaperSize',[1400 500]); hold on
subplot(1,3,1), hold on
xlabel('Time (s)')
ylabel('Firing rate (norm)')
c = parula(length(paradigm_names));
for i = 1:length(paradigm_names)
	for j = 1:max(combined_data.experiment)
		subplot(1,3,j), hold on
		plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
		plot_these = intersect(plot_these, find(combined_data.experiment==j));

		if ~isempty(plot_these)
			this_resp=mean2(combined_data.fA(:,plot_these));
			time = dt*(1:length(this_resp));
			time = time(a:z);
			this_resp = this_resp(a:z);
			this_resp = this_resp/max(this_resp);
			plot(time,this_resp,'Color',c(i,:))
		end
	end
end

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

% Does a fractional differentiation model explain this response kinetics? How does the order of the fractional differnetiator vary? 

% a = floor(5/dt);
% z = floor(12/dt);
% if redo
% 	clear fdmodel
% 	for i = 1:length(paradigm_names)
% 		fdmodel(i).d = [];
% 		fdmodel(i).alpha = [];
% 		fdmodel(i).fp = [];
% 		plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
% 		this_resp=mean2(combined_data.fA(:,plot_these));
% 		this_pid=mean2(combined_data.PID(plot_these,:));
% 		time = dt*(1:length(this_resp));
% 		time = time(a:z);
% 		this_resp = this_resp(a:z);
% 		this_pid = this_pid(a:z);
% 		[fdmodel(i).alpha,fdmodel(i).d,fdmodel(i).fp] =  FitFractionalDModel(this_pid, this_resp,100);
% 	end
% end


% figure('outerposition',[0 0 1200 800],'PaperUnits','points','PaperSize',[1200 800]); hold on
% for i = 1:length(paradigm_names)
% 	subplot(2,4,i), hold on
% 	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
% 	this_resp=mean2(combined_data.fA(:,plot_these));
% 	time = dt*(1:length(this_resp));
% 	this_resp = this_resp(a:z);
% 	time = time(a:z);
% 	plot(time,this_resp,'k')
% 	plot(time,fdmodel(i).fp,'r')
% 	title(strcat('alpha = ',oval(fdmodel(i).alpha,2)),'interpreter','tex')
% 	xlabel('Time (s)')
% 	ylabel('Firing Rate (Hz)')
% end


% PrettyFig;
% if being_published
% 	snapnow
% 	delete(gcf)
% end

%%
% It looks like in an effort to get the adaptation ratio right, the fractional diff. models are messing up the timescales. 


%% Very Long Steps
% In this section, we look at responses of ORNs to very long steps (2min+). The following figure shows the raw traces: 

load('very_long_steps.mat')

dt = 3e-3;
a = floor(2/dt);
z = floor(152/dt);
z1 = floor(151/dt);
max_pid = [];
max_f = [];
end_pid = [];
end_f = [];

figure('outerposition',[0 0 1500 1000],'PaperUnits','points','PaperSize',[1500 1000]); hold on
subplot(2,1,1), hold on
c = parula(length(paradigm_names));
for i = 1:length(paradigm_names)
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	this_pid=mean2(combined_data.PID(plot_these,:));
	time = dt*(1:length(this_pid));
	end_pid = [end_pid mean(this_pid(z1:z))];
	time = time(a:z);
	this_pid = this_pid(a:z);
	max_pid = [max_pid max(this_pid)];
	plot(time,this_pid,'Color',c(i,:))
end
xlabel('Time (s)')
ylabel('Stimulus (V)')
set(gca,'YScale','log')

subplot(2,1,2), hold on
for i = 1:length(paradigm_names)
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	this_resp=mean2(combined_data.fA(:,plot_these));
	time = dt*(1:length(this_resp));
	end_f = [end_f mean(this_resp(z1:z))];
	time = time(a:z);
	this_resp = this_resp(a:z);
	max_f = [max_f max(this_resp)];
	plot(time,this_resp,'Color',c(i,:))
end
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')


PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%%
% The following figure shows how the response and stimulus are related, for the peak, and for the terminal 1 second. 

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(max_pid,max_f,'r')
plot(end_pid,end_f,'k')
legend({'Peak','Terminal'},'Location','NorthWest')
xlabel('Stimulus (V)')
ylabel('Firing Rate (Hz)')

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end


%%
% What does the initial transient look like? 


a = floor(4/dt);
z = floor(12/dt);


figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
c = parula(length(paradigm_names));
for i = 1:length(paradigm_names)
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	this_pid=mean2(combined_data.PID(plot_these,:));
	time = dt*(1:length(this_pid));
	time = time(a:z);
	this_pid = this_pid(a:z);
	plot(time,this_pid,'Color',c(i,:))
end
xlabel('Time (s)')
ylabel('Stimulus (V)')
set(gca,'YScale','log')

L ={};
lh=[];
subplot(1,2,2), hold on
for i = 1:length(paradigm_names)
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	this_resp=mean2(combined_data.fA(:,plot_these));
	time = dt*(1:length(this_resp));
	time = time(a:z);
	this_resp = this_resp(a:z);
	if i > 2
		lh = [lh plot(time,this_resp,'Color',c(i,:))];
		[~,loc]=max(this_resp);
		t = time(loc:end);
		t = t -min(t);
		ft=fittype('a*exp(-x./b1)+ a*exp(-x./b2) +c');
		fo = fitoptions(ft);
		fo.Robust = 'on';
		pf = max(this_resp) -  this_resp(end); pf = pf/2;
		fo.Lower = [pf-1e-3 1e-2 1e-2 1];
		fo.Upper = [pf+1e-3 2 20 max(this_resp)];
		fo.StartPoint = [max(this_resp) .1 1 min(this_resp)];
		ff = fit(t(:),this_resp(loc:end),ft,fo);
		plot(time(loc:end),ff(t),'r')
		L = [L strcat('\tau=',oval(ff.b1,2),',',oval(ff.b2,2),'s'); ];
	else
		plot(time,this_resp,'Color',c(i,:))
	end
end
xlabel('Time (s)')
legend(lh,L)
ylabel('Firing Rate (Hz)')

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%% When does the response stop adapting?
% When does the ORN response reach steady state? To determine this, we chunk the data into 5-second blocks, and fit lines and extract slopes from each chunk, and plot them as a function of time. 



% plot_data is indexed by where we start
all_start = [5:5:145];
all_end = all_start+5;

clear plot_data
for i = 3:length(paradigm_names)
	plot_data(i).resp_slope = [];
	plot_data(i).resp_slope_err = [];
	plot_data(i).resp_mean = [];
	plot_data(i).resp_mean_err = [];

	for j = 1:length(all_start)
		a = floor(all_start(j)/dt);
		z = floor(all_end(j)/dt);
		n = sqrt(z-a);

		plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
		these_resp=mean2(combined_data.fA(:,plot_these));

		cropped_resp = these_resp(a:z);
		t = dt*(1:length(cropped_resp));

		plot_data(i).resp_mean = 		[plot_data(i).resp_mean mean(cropped_resp)];
		plot_data(i).resp_mean_err = 	[plot_data(i).resp_mean_err std(cropped_resp)/n];

		ff = fit(t(:),cropped_resp(:),'poly1');
		gof=confint(ff);
		gof=gof(:,1);

		plot_data(i).resp_slope = [plot_data(i).resp_slope ff.p1];
		plot_data(i).resp_slope_err = [plot_data(i).resp_slope_err gof(2)-gof(1)];

	end
end

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
c = parula(length(paradigm_names));
for i = 3:length(plot_data)
	errorbar(all_start+2.5,plot_data(i).resp_slope,plot_data(i).resp_slope_err,'Color',c(i,:))

end
plot([0 all_start(end)],[0 0],'--')
xlabel('Time (s)')
ylabel('Slope (Hz/s)')
L = paradigm_names(3:4);
for i = 1:length(L)
	L{i} = strrep(L{i},'_','-'); 
end
legend(L)

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%% What is the functional form of the decay in response?
% So far, we have been fitting an exponential to the response in order to extract a timescale. But we didn't claim that it was actually an exponential. Now that we have a very long trajectory, we will fit various functional forms to this and see which one fits best. 

warning off
lh=[];
a = floor(2/dt);
z = floor(150/dt);
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
for i = 1:2

	subplot(1,2,i), hold on
	plot_these=find(strcmp(paradigm_names{2+i}, combined_data.paradigm));
	this_resp=mean2(combined_data.fA(:,plot_these));
	this_resp = this_resp(a:z);
	time = dt*(1:length(this_resp));
	[~,loc]=max(this_resp);
	t = time(loc:end);
	t = t -min(t) + dt;
	plot(t,this_resp(loc:end),'k')

	ft=fittype('a*exp(-x./b)+c');
	options = fitoptions(ft);
	options.StartPoint = [max(this_resp) .4 10];
	ff = fit(t(:),this_resp(loc:end),ft,options);
	lh(1)=plot(t,ff(t),'r');
	allfits(i,1).ff = ff;

	ft=fittype('a*exp(-x./b)+ a2*exp(-x./b2) + c');
	options = fitoptions(ft);
	options.StartPoint = [max(this_resp) max(this_resp) .4 .4 10];
	options.Upper = [max(this_resp) max(this_resp) 10 10 max(this_resp)];
	options.Lower = [10 10 dt dt 10];
	ff = fit(t(:),this_resp(loc:end),ft,options);
	lh(2)=plot(t,ff(t),'b');
	allfits(i,2).ff = ff;

	ft=fittype('a*exp(-x./b)+ a2*exp(-x./b2) + a3*exp(-x./b3)  + c');
	options = fitoptions(ft);
	options.StartPoint = [max(this_resp) max(this_resp)  max(this_resp) .4 .4 .4 10];
	options.Upper = [max(this_resp) max(this_resp) max(this_resp) 10 10 10 max(this_resp)];
	options.Lower = [10 10 10 dt dt dt 10];
	ff = fit(t(:),this_resp(loc:end),ft,options);
	lh(3)=plot(t,ff(t),'m');
	allfits(i,3).ff = ff;

	ft=fittype('a*(x.^b) + c');
	options = fitoptions(ft);
	options.StartPoint = [max(this_resp) -1 10];
	options.Upper = [max(this_resp) -.1 max(this_resp)];
	options.Lower = [10 -100 10];
	ff = fit(t(:),this_resp(loc:end),ft,options);
	lh(4)=plot(t,ff(t),'g');
	allfits(i,4).ff = ff;

	set(gca,'YLim',[min(this_resp(loc:end))-2 this_resp(loc)*2])
	set(gca,'XScale','log');
	set(gca,'YScale','log')
	L = {'exp1','exp2','exp3','power law'};
	legend(lh,L,'Location','SouthWest')
	xlabel('Time (s)')
	ylabel('Firing rate (Hz)')
	title(strrep(paradigm_names{i+2},'_','-')) 
end

warning on

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%%
% The exponents of the power law fits in the two cases are:

disp([allfits(1,4).ff.b allfits(2,4).ff.b])



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


