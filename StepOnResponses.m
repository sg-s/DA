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
% % What do the responses of the ORN to the step ON of the odor look like? The following figure shows the stimulus on the left and the response on the right, for the various heights of the step ON. 

load('MeanShiftedGaussians.mat')


a = floor(4/3e-3);
z = floor(12/3e-3);
z1 = floor(11/3e-3);
max_pid = [];
max_f = [];
end_pid = [];
end_f = [];

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1), hold on
c = parula(length(paradigm_names));
for i = 1:length(paradigm_names)
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	this_pid=mean2(combined_data.PID(plot_these,:));
	time = 3e-3*(1:length(this_pid));
	end_pid = [end_pid mean(this_pid(z1:z))];
	time = time(a:z);
	this_pid = this_pid(a:z);
	max_pid = [max_pid max(this_pid)];
	plot(time,this_pid,'Color',c(i,:))
end
xlabel('Time (s)')
ylabel('Stimulus (V)')


subplot(1,3,2), hold on
for i = 1:length(paradigm_names)
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	this_resp=mean2(combined_data.fA(:,plot_these));
	time = 3e-3*(1:length(this_resp));
	end_f = [end_f mean(this_resp(z1:z))];
	time = time(a:z);
	this_resp = this_resp(a:z);
	plot(time,this_resp,'Color',c(i,:))
	max_f = [max_f max(this_resp)];

end
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')

subplot(1,3,3), hold on
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
% What are the timescales of relaxation of firing rates in each of these cases? In the following figure, we fit an exponential to the firing rates from peak onwards, not to claim that the decay is exponential, but to extract a timescale for each case. 

figure('outerposition',[0 0 1200 800],'PaperUnits','points','PaperSize',[1200 800]); hold on
for i = 1:length(paradigm_names)
	subplot(2,3,i), hold on
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	this_resp=mean2(combined_data.fA(:,plot_these));
	time = 3e-3*(1:length(this_resp));
	time = time(a:z);
	this_resp = this_resp(a:z);
	[~,loc]=max(this_resp);
	t = time(loc:end);
	t = t -min(t);
	ft=fittype('a*exp(-x./b)+c');
	warning off
	ff = fit(t(:),this_resp(loc:end),ft);
	warning on
	plot(time,this_resp,'k')
	plot(time(loc:end),ff(t),'r')
	title(strcat('tau=',oval(ff.b,2),'s'),'interpreter','tex')
	xlabel('Time (s)')
	ylabel('Firing Rate (Hz)')
end

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%%
% Does a fractional differentiation model explain this response kinetics? How does the order of the fractional differnetiator vary? 

a = floor(5/3e-3);
z = floor(12/3e-3);
if redo
	clear fdmodel
	for i = 1:length(paradigm_names)
		fdmodel(i).d = [];
		fdmodel(i).alpha = [];
		fdmodel(i).fp = [];
		plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
		this_resp=mean2(combined_data.fA(:,plot_these));
		this_pid=mean2(combined_data.PID(plot_these,:));
		time = 3e-3*(1:length(this_resp));
		time = time(a:z);
		this_resp = this_resp(a:z);
		this_pid = this_pid(a:z);
		[fdmodel(i).alpha,fdmodel(i).d,fdmodel(i).fp] =  FitFractionalDModel(this_pid, this_resp,100);
	end
end


figure('outerposition',[0 0 1200 800],'PaperUnits','points','PaperSize',[1200 800]); hold on
for i = 1:length(paradigm_names)
	subplot(2,3,i), hold on
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	this_resp=mean2(combined_data.fA(:,plot_these));
	time = 3e-3*(1:length(this_resp));
	this_resp = this_resp(a:z);
	time = time(a:z);
	plot(time,this_resp,'k')
	plot(time,fdmodel(i).fp,'r')
	title(strcat('alpha = ',oval(fdmodel(i).alpha,2)),'interpreter','tex')
	xlabel('Time (s)')
	ylabel('Firing Rate (Hz)')
end


PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%%
% It looks like in an effort to get the adaptation ratio right, the fractional diff. models are messing up the timescales. 


%% Very Long Steps
% In this section, we look at responses of ORNs to very long steps (2min+). The following figure shows the raw traces: 

load('very_long_steps.mat')


a = floor(2/3e-3);
z = floor(152/3e-3);
z1 = floor(151/3e-3);
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
	time = 3e-3*(1:length(this_pid));
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
	time = 3e-3*(1:length(this_resp));
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


a = floor(4/3e-3);
z = floor(12/3e-3);


figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
c = parula(length(paradigm_names));
for i = 1:length(paradigm_names)
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	this_pid=mean2(combined_data.PID(plot_these,:));
	time = 3e-3*(1:length(this_pid));
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
	time = 3e-3*(1:length(this_resp));
	time = time(a:z);
	this_resp = this_resp(a:z);
	if i > 2
		lh = [lh plot(time,this_resp,'Color',c(i,:))];
		[~,loc]=max(this_resp);
		t = time(loc:end);
		t = t -min(t);
		ft=fittype('a*exp(-x./b)+c');
		options = fitoptions(ft);
		options.StartPoint = [max(this_resp) .4 10];
		ff = fit(t(:),this_resp(loc:end),ft,options);
		plot(time(loc:end),ff(t),'r')
		L = [L strcat('\tau=',oval(ff.b,2),'s')];
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
		a = floor(all_start(j)/3e-3);
		z = floor(all_end(j)/3e-3);
		n = sqrt(z-a);

		plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
		these_resp=mean2(combined_data.fA(:,plot_these));

		cropped_resp = these_resp(a:z);
		t = 3e-3*(1:length(cropped_resp));

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
	errorbar(all_start+2.5,plot_data(i).resp_slope,plot_data(i).resp_slope_err)

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
a = floor(2/3e-3);
z = floor(150/3e-3);
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
for i = 1:2

	subplot(1,2,i), hold on
	plot_these=find(strcmp(paradigm_names{2+i}, combined_data.paradigm));
	this_resp=mean2(combined_data.fA(:,plot_these));
	this_resp = this_resp(a:z);
	time = 3e-3*(1:length(this_resp));
	[~,loc]=max(this_resp);
	t = time(loc:end);
	t = t -min(t) + 3e-3;
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
	options.Lower = [10 10 3e-3 3e-3 10];
	ff = fit(t(:),this_resp(loc:end),ft,options);
	lh(2)=plot(t,ff(t),'b');
	allfits(i,2).ff = ff;

	ft=fittype('a*exp(-x./b)+ a2*exp(-x./b2) + a3*exp(-x./b3)  + c');
	options = fitoptions(ft);
	options.StartPoint = [max(this_resp) max(this_resp)  max(this_resp) .4 .4 .4 10];
	options.Upper = [max(this_resp) max(this_resp) max(this_resp) 10 10 10 max(this_resp)];
	options.Lower = [10 10 10 3e-3 3e-3 3e-3 10];
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


