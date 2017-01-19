
pHeader;


%% Estimating the timescale of gain control
% In this document, I attempt to determine the timescale of gain control from data where we present pulses of odorant on top of a background, at various time-points since background onset. 

% first, gather the data
root = [getPath(dataManager,'b1840fa24e8f3070315f6a7e83358d4e') oss];
allfiles = dir([root '*.kontroller']);
allfiles(cellfun(@(x) strcmp(x(1),'.'), {allfiles.name})) = [];

S = [];
fA = [];
LFP = [];
fly = [];
sensillum = [];
ab = [];
pulse_time = [];

for i = 1:length(allfiles)
	clear data spikes ControlParadigm
	load([root allfiles(i).name],'-mat')

	this_fly = str2double(allfiles(i).name(strfind(allfiles(i).name,'_F')+2));
	this_ab = str2double(allfiles(i).name(strfind(allfiles(i).name,'_ab')+3));
	this_sensillum = str2double(allfiles(i).name(strfind(allfiles(i).name,'_S')+2));

	this_S = vertcat(data.PID)'; this_S = this_S(1:10:end,:);
	this_LFP = vertcat(data.voltage)'; this_LFP = this_LFP(1:10:end,:);

	this_fA = spiketimes2f(vertcat(spikes.A),1e-4*(1:length(vertcat(spikes.A))),1e-3,3e-2);

	this_fly = this_fly*ones(size(this_S,2),1);
	this_sensillum = this_sensillum*ones(size(this_S,2),1);
	this_ab = this_ab*ones(size(this_S,2),1);

	% figure out the time of the pulse 
	n_trials_per_paradigm = cellfun(@width,{data.PID});
	this_pulse_time = [];
	for j = 1:length(ControlParadigm)
		tpt = find(ControlParadigm(j).Outputs(5,:),1,'first')*1e-4 - find(ControlParadigm(j).Outputs(4,:),1,'first')*1e-4;
		this_pulse_time = [this_pulse_time; zeros(n_trials_per_paradigm(j),1)+tpt];
	end

	% consolidate
	fly = [fly; this_fly];
	sensillum = [sensillum; this_sensillum];
	pulse_time = [pulse_time; this_pulse_time];
	S = [S this_S];
	fA = [fA this_fA];
	LFP = [LFP this_LFP];
	ab = [ab; this_ab];
end

% remove some crappy trials
rm_this = isnan(sum(LFP)) == 1;
LFP(:,rm_this) = [];
S(:,rm_this) = [];
fA(:,rm_this) = [];
pulse_time(rm_this) = [];
sensillum(rm_this) = [];
fly(rm_this) = [];
ab(rm_this) = [];

% remove baselines from LFP and stimulus
for i = 1:length(sensillum)
	LFP(:,i) = LFP(:,i) - mean(LFP(1:1e3,i));
	S(:,i) = S(:,i) - min(S(1:end,i));
end

time = 1e-3*(1:length(S));

%%
% First, I plot the raw stimulus and the LFP responses, so we can get a sense of the data. In the following figure, each column is data from a single neuron (the identity of which is indicated in the title). The top row shows the stimulus, the middle row shows the LFP, and the bottom row shows the firing rate. Many firing rate traces are missing because I didn't sort of all the data. Note that the stimulus seems quite reliable, with no systematic trend. Also note that the LFP responses to the stimulus pulse seem identical -- independent of the time since the step turns on.

%%
% This is quite unexpected. Note too that we tried different background heights, foreground pulse heights. Also, because we did it on ab3A and ab2A, and because ab2A is much more sensitive to this odorant than ab3A, we have also done this experiment effectively at two very different concentration regimes: close to saturation (for ab2A), and in the midpoint of the sensitivity range (for ab3A). 


all_sensillum = unique(sensillum);

figure('outerposition',[0 0 1411 704],'PaperUnits','points','PaperSize',[1411 704]); hold on
for si = 1:length(all_sensillum)
	ts = all_sensillum(si);

	% plot stimulus
	subplot(3,length(all_sensillum),si); hold on
	all_pulse_times = unique(pulse_time);
	c = parula(length(all_pulse_times)+1);
	for i = 1:length(all_pulse_times)
		y = mean(S(:,pulse_time == all_pulse_times(i) & sensillum == ts),2);
		plot(time(1:50:end),y(1:50:end),'Color',c(i,:))
	end
	if si == 1
		ylabel('Stimulus (V)')
	end
	set(gca,'XLim',[0 10])
	title(['ab' oval(mean(ab(sensillum == ts)))])

	% plot LFP
	subplot(3,length(all_sensillum),length(all_sensillum)+si); hold on
	all_pulse_times = unique(pulse_time);
	c = parula(length(all_pulse_times)+1);
	for i = 1:length(all_pulse_times)
		y = mean(LFP(:,pulse_time == all_pulse_times(i) & sensillum == ts),2);
		plot(time(1:50:end),y(1:50:end),'Color',c(i,:))
	end
	if si == 1
		ylabel('\Delta LFP (V)')
	end
	set(gca,'XLim',[0 10])

	% plot firing rate
	subplot(3,length(all_sensillum),2*length(all_sensillum)+si); hold on
	all_pulse_times = unique(pulse_time);
	c = parula(length(all_pulse_times)+1);
	for i = 1:length(all_pulse_times)
		y = fA(:,pulse_time == all_pulse_times(i) & sensillum == ts);
		y(:,sum(y)==0) = [];
		y = mean(y,2);
		if ~isempty(y)
			plot(time(1:50:end),y(1:50:end),'Color',c(i,:))
		end
	end
	xlabel('Time (s)')
	if si  == 1
		ylabel('Firing rate (Hz)')
	end
	set(gca,'XLim',[0 10])

end

prettyFig();

if being_published
	snapnow
	delete(gcf)
end


%%
% Now, I quantify how these responses vary as a function of time since background odour onset. In the following figure, each column is data from a single neuron, the identity of which is indicated in the title. The first row shows the peak stimulus, the second row the peak LFP response. The third row shows the change in the LFP due to the stimulus, and the fourth row is the effective gain of the LFP (change in LFP response/change in stimulus). For each plot, I also calculate the Spearman correlation coefficient ($\rho$). 

%%
% Naively, I would expect that the gain decreases with time -- as the neuron adapts to the background pulse, the gain becomes smaller and smaller. However, I see an *increase* in gain in 4/5 neurons (the exception being the fourth column). In all cases, the change in the gain seems small and inconsequential. Does this mean the neuron is not adapting at all? 

% compute peak_LFP, etc
peak_LFP = NaN*pulse_time;
delta_LFP = NaN*pulse_time;
delta_S = NaN*pulse_time;

for i = 1:length(pulse_time)
	a = (1 + pulse_time(i))*1e3;
	z = a + 1e3;
	peak_LFP(i) = min(LFP(a:z,i));
	delta_LFP(i) = peak_LFP(i) - mean(LFP(a:a+100,i));
	delta_S(i) = max(S(:,i)) - mean(S(a:a+100,i));
end

figure('outerposition',[0 0 1411 901],'PaperUnits','points','PaperSize',[1411 901]); hold on

marker_colour = .2*ones(3,1);

for si = 1:length(all_sensillum)
	ts = all_sensillum(si);
	x = pulse_time(sensillum == ts);


	subplot(4,length(all_sensillum),si); hold on
	y = max(S(:,sensillum == ts));
	plot(x,y,'.','MarkerSize',24,'Color',marker_colour)
	set(gca,'XLim',[0 5],'YLim',[0 max(y)*2],'XTick',0:5)
	legend({['\rho = ' oval(spear(x(:),y(:)))]},'Location','northeast')
	if si == 1
		ylabel('S_{max}')
	end
	title(['ab' oval(mean(ab(sensillum == ts)))])

	subplot(4,length(all_sensillum),length(all_sensillum)+si); hold on
	y = peak_LFP(sensillum == ts);
	plot(x,y,'.','MarkerSize',24,'Color',marker_colour)
	if si == 1
		ylabel('LFP_{peak}')
	end
	set(gca,'XLim',[0 5],'YLim',[-3 0],'XTick',0:5)
	legend({['\rho = ' oval(spear(x(:),y(:)))]},'Location','southeast')


	subplot(4,length(all_sensillum),2*length(all_sensillum)+si); hold on
	y = -delta_LFP(sensillum == ts);
	plot(x,y,'.','MarkerSize',24,'Color',marker_colour)
	if si == 1
		ylabel('abs(\DeltaLFP)')
	end
	legend({['\rho = ' oval(spear(x(:),y(:)))]},'Location','northeast')
	set(gca,'XLim',[0 5],'YLim',[0 1],'XTick',0:5)

	subplot(4,length(all_sensillum),3*length(all_sensillum)+si); hold on
	y = -delta_LFP(sensillum == ts)./delta_S(sensillum == ts);
	plot(x,y,'.','MarkerSize',24,'Color',marker_colour)
	if si == 1
		ylabel('abs(\DeltaLFP)/\DeltaS')
	end
	xlabel('Time since odour onset (s)')
	legend({['\rho = ' oval(spear(x(:),y(:)))]},'Location','northeast')
	set(gca,'XLim',[0 5],'YLim',[0 3*max(y)],'XTick',0:5)


end

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%% Reproducing LFP responses precisely
% In this section, I attempt to fit various binding models to the LFP traces, to try to understand what is going on more quantiatively, and also to see if this data constrains the space of models in any way.

%%
% In the following figure, I fit various models to LFP responses from the ab2 sensillum. I fit only two traces from a single neuron here. Each panel shows the ODEs of the model being fit, the responses of the sensillum in black, and the best-fit model predictions in red. In addition, I show some time-scale parameters where possible. The first panel is the fit from a non-adapting model; all other panels show fits from models that adapt in various ways. Note that in every case, it appears as though the kinetics of the response tot he pulse in the LFP data appear to be slower than the kinetics of the model response. Also note that no model appears to get the heights of responses to both pulses correctly. 

clear data
all_pulse_times = unique(pulse_time);
for i = 1:length(all_pulse_times)
	tpt = all_pulse_times(i);
	j = find(pulse_time == tpt & sensillum == 3,1,'first');
	data(i).stimulus = S(:,j);
	data(i).response = -LFP(:,j);
	data(i).response(1:1e3) = NaN;
end



figure('outerposition',[0 0 1300 900],'PaperUnits','points','PaperSize',[1300 900]); hold on

% no adaptation model (LFPmodelv1)
clear p
p.      k1 = 1.5886e+03;
p.      k2 = 42.9761;
p. R_scale = 2.4621;
p.R_offset = -0.1412;

subplot(2,3,1); hold on
for i = [1 9]
	plot([data(i).response],'k')
	R = LFPmodelv1(data(i).stimulus,p);
	plot(R,'r')
end
xlabel('Time (ms)')
set(gca,'XLim',[1e3 6e3])
ylabel('\DeltaLFP (mV)')
text(2000,2.4,'$\dot{a}=k_{+}(1-a)S-k_{-}a$','Interpreter','latex','FontSize',20)
 

% adapting model, k_D changes with S, k2 fixed
clear p
p.adap_tau = 4.6875;
p.      k2 = 6.2150;
p.   k_min = 0.0131;
p. R_scale = 2.6319;
p.R_offset = -0.1221;

subplot(2,3,2); hold on
for i = [1 9]
	plot([data(i).response],'k')
	R = LFPmodelv2(data(i).stimulus,p);
	plot(R,'r')
end
xlabel('Time (ms)')
set(gca,'XLim',[1e3 6e3])
text(2000,1,'$\tau_{A}\dot{k_{D}}=S-k_{D}$','Interpreter','latex','FontSize',20)
text(2000,.5,'$k_{+}=k_{-}/k_{D}$','Interpreter','latex','FontSize',20)
text(2000,0,['$\tau_A$ = ' oval(p.adap_tau) 's'],'Interpreter','latex','FontSize',20)
ylabel('\DeltaLFP (mV)')

% adapting model, k_D changes with S, k1 fixed 
clear p
p.adap_tau = 5.8672;
p.   k_min = 0.0116;
p. R_scale = 2.5303;
p.R_offset = -0.1192;
p.      k1 = 331;
subplot(2,3,3); hold on
for i = [1 9]
	plot([data(i).response],'k')
	R = LFPmodelv3(data(i).stimulus,p);
	plot(R,'r')
end
xlabel('Time (ms)')
set(gca,'XLim',[1e3 6e3])
text(2000,1,'$\tau_{A}\dot{k_{D}}=S-k_{D}$','Interpreter','latex','FontSize',20)
text(2000,.5,'$k_{-}=k_{+}k_{D}$','Interpreter','latex','FontSize',20)
text(2000,0,['$\tau_A$ = ' oval(p.adap_tau) 's'],'Interpreter','latex','FontSize',20)
ylabel('\DeltaLFP (mV)')

% LFP model v4
% adapting model, k_D changes like in chemotaxis, k1 fixed
clear p
p.adap_tau = 1.3750;
p.   k_min = 0.0108;
p. R_scale = 2.4210;
p.R_offset = -0.1875;
p.      k1 = 539;
subplot(2,3,4), hold on
for i = [1 9]
	plot([data(i).response],'k')
	R = LFPmodelv4(data(i).stimulus,p);
	plot(R,'r')
end
xlabel('Time (ms)')
set(gca,'XLim',[1e3 6e3])
text(2000,1,'$\tau_{A}\dot{k_{D}}=k_{D}(a-1/2)$','Interpreter','latex','FontSize',20)
text(2000,.5,'$k_{-}=k_{+}k_{D}$','Interpreter','latex','FontSize',20)
text(2000,0,['$\tau_A$ = ' oval(p.adap_tau) 's'],'Interpreter','latex','FontSize',20)
ylabel('\DeltaLFP (mV)')

% LFP model v5
% adapting model, K_D changes as in chemotaxis, k2 fixed
clear p
p.adap_tau=  1.0547;
p.   k_min=  0.0089;
p. R_scale=  2.4835;
p.R_offset=  -0.2036;
p.      k2=  5.5937;
subplot(2,3,5), hold on
for i = [1 9]
	plot([data(i).response],'k')
	R = LFPmodelv5(data(i).stimulus,p);
	plot(R,'r')
end
xlabel('Time (ms)')
set(gca,'XLim',[1e3 6e3])
text(2000,1,'$\tau_{A}\dot{k_{D}}=k_{D}(a-1/2)$','Interpreter','latex','FontSize',20)
text(2000,.5,'$k_{+}=k_{-}/k_{D}$','Interpreter','latex','FontSize',20)
text(2000,0,['$\tau_A$ = ' oval(p.adap_tau) 's'],'Interpreter','latex','FontSize',20)
ylabel('\DeltaLFP (mV)')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%% Non-saturated responses
% One concern in the previous dataset is that the probe pulses can be saturating, and this could be the reason why we don't see any change. So we repeated the experiment with smaller pulses, and even went into a regime where the probe pulse was so small w.r.t to the background that we couldn't see it (but the neuron could). 

%%
% The following figure shows the responses of three neurons to pulses at various times. There are four columns because we repeated the experiment with two different values of the background stimulus. Note that in every case, the LFP and firing rate responses to the probe pulse are almost identical, independent of the location of the probe pulse relative to the step on in the stimulus background. Note too that in all cases, the neuron is far from saturation, yet the responses are extremely stereotyped. 

clearvars -except being_published
% first, gather the data
root = [getPath(dataManager,'fdd9f9238e58fceac600d9b45f6122da') oss];
allfiles = dir([root '*.kontroller']);
allfiles(cellfun(@(x) strcmp(x(1),'.'), {allfiles.name})) = [];

S = [];
fA = [];
LFP = [];
fly = [];
sensillum = [];
ab = [];
pulse_time = [];

for i = 1:length(allfiles)
	clear data spikes ControlParadigm
	load([root allfiles(i).name],'-mat')

	this_fly = str2double(allfiles(i).name(strfind(allfiles(i).name,'_F')+2));
	this_ab = str2double(allfiles(i).name(strfind(allfiles(i).name,'_ab')+3));
	this_sensillum = str2double(allfiles(i).name(strfind(allfiles(i).name,'_S')+2));

	this_S = vertcat(data.PID)'; this_S = this_S(1:10:end,:);
	this_LFP = vertcat(data.voltage)'; this_LFP = this_LFP(1:10:end,:);

	this_fA = spiketimes2f(vertcat(spikes.A),1e-4*(1:length(vertcat(spikes.A))),1e-3,3e-2);

	this_fly = this_fly*ones(size(this_S,2),1);
	this_sensillum = this_sensillum*ones(size(this_S,2),1);
	this_ab = this_ab*ones(size(this_S,2),1);

	% figure out the time of the pulse 
	n_trials_per_paradigm = cellfun(@width,{data.PID});
	this_pulse_time = [];
	for j = 1:length(ControlParadigm)
		tpt = find(ControlParadigm(j).Outputs(5,:),1,'first')*1e-4 - find(ControlParadigm(j).Outputs(4,:),1,'first')*1e-4;
		this_pulse_time = [this_pulse_time; zeros(n_trials_per_paradigm(j),1)+tpt];
	end

	% consolidate
	fly = [fly; this_fly];
	sensillum = [sensillum; this_sensillum];
	pulse_time = [pulse_time; this_pulse_time];
	S = [S this_S];
	fA = [fA this_fA];
	LFP = [LFP this_LFP];
	ab = [ab; this_ab];
end

% remove some crappy trials
rm_this = isnan(sum(LFP)) == 1;
LFP(:,rm_this) = [];
S(:,rm_this) = [];
fA(:,rm_this) = [];
pulse_time(rm_this) = [];
sensillum(rm_this) = [];
fly(rm_this) = [];
ab(rm_this) = [];

% remove baselines from LFP and stimulus
for i = 1:length(sensillum)
	LFP(:,i) = LFP(:,i) - mean(LFP(1:1e3,i));
	S(:,i) = S(:,i) - min(S(1:end,i));
end

time = 1e-3*(1:length(S));

all_sensillum = unique(sensillum);

figure('outerposition',[0 0 1411 704],'PaperUnits','points','PaperSize',[1411 704]); hold on
for si = 1:length(all_sensillum)
	ts = all_sensillum(si);

	% plot stimulus
	subplot(3,length(all_sensillum),si); hold on
	all_pulse_times = unique(pulse_time);
	c = parula(length(all_pulse_times)+1);
	for i = 1:length(all_pulse_times)
		y = mean(S(:,pulse_time == all_pulse_times(i) & sensillum == ts),2);
		plot(time(1:50:end),y(1:50:end),'Color',c(i,:))
	end
	if si == 1
		ylabel('Stimulus (V)')
	end
	set(gca,'XLim',[0 10])
	title(['ab' oval(mean(ab(sensillum == ts)))])

	% plot LFP
	subplot(3,length(all_sensillum),length(all_sensillum)+si); hold on
	all_pulse_times = unique(pulse_time);
	c = parula(length(all_pulse_times)+1);
	for i = 1:length(all_pulse_times)
		y = mean(LFP(:,pulse_time == all_pulse_times(i) & sensillum == ts),2);
		plot(time(1:50:end),y(1:50:end),'Color',c(i,:))
	end
	if si == 1
		ylabel('\Delta LFP (V)')
	end
	set(gca,'XLim',[0 10])

	% plot firing rate
	subplot(3,length(all_sensillum),2*length(all_sensillum)+si); hold on
	all_pulse_times = unique(pulse_time);
	c = parula(length(all_pulse_times)+1);
	for i = 1:length(all_pulse_times)
		y = fA(:,pulse_time == all_pulse_times(i) & sensillum == ts);
		y(:,sum(y)==0) = [];
		y = mean(y,2);
		if ~isempty(y)
			plot(time(1:50:end),y(1:50:end),'Color',c(i,:))
		end
	end
	xlabel('Time (s)')
	if si  == 1
		ylabel('Firing rate (Hz)')
	end
	set(gca,'XLim',[0 10])

end

prettyFig();

if being_published
	snapnow
	delete(gcf)
end



%% Estimating timescale of gain control from Gaussian data
% Can I estimate timescales of gain control from the Gaussian data? In the following figure, I fit exponentials to the firing rates from the Gaussian data, and subtract the responses from these exponential fits. I then estimate the absolute deivation from these expoenentials, and fit another exponential to these absolute deviations, which is my proxy for the timescale of gain control. 

clearvars -except being_published

% define what we want to work on
data_hashes = {'93ba5d68174e3df9f462a1fc48c581da','bcd4cf4fe12817d084a2b06f981161ee','cd6753c0e4cf02895cd5e2c5cb58aa1a','3ea08ccfa892c6545d74bbdaaa6cbee1','f11c4a5792d0c9fec7c40fd6aa2fce40'};
odour_names = {'ethyl-acetate','1-pentanol','1-pentanol','2-butanone','isoamyl-acetate'};
orn_names = {'ab3A','ab3A','ab2A','ab2A','pb1A'};

if exist('.cache/tau_gain_MSG.mat','file') == 2
	load('.cache/tau_gain_MSG.mat','allmetrics')
else
	allmetrics = struct([]);

	for di = 1:length(data_hashes)
		clear MSGdata
		MSGdata = consolidateData2(getPath(dataManager,data_hashes{di}));

		% make sure stimulus is always positive
		for i = 1:size(MSGdata.PID,2)
			MSGdata.PID(:,i) = MSGdata.PID(:,i) - min(MSGdata.PID(:,i));
		end


		% look at timescale of response, fluctuations
		allmetrics(di).fA_tau = NaN*MSGdata.fly;
		allmetrics(di).fA_gain_tau = NaN*MSGdata.fly;
		allmetrics(di).fA_gain_tau_err = NaN*MSGdata.fly;
		allmetrics(di).S_tau = NaN*MSGdata.fly;
		allmetrics(di).orn = MSGdata.orn;
		ft = fittype('A*exp(-x/tau) + B');
		ft_stim = fittype('A*(1- exp(-x/tau))');
		x = 1:50e3;
		for i = 1:width(MSGdata.PID)
			textbar(i,width(MSGdata.PID))
			R = MSGdata.fA(5e3+1:55e3,i);
			if max(R)>0
				ff = fit(x(:),R,ft,'StartPoint',[max(R) 0 1e3],'Lower',[1 0 10]);
				allmetrics(di).fA_tau(i) = ff.tau;
				y = abs(R - ff(x));
				ff = fit(x(:),y,ft,'StartPoint',[max(y) 0 1e2],'Lower',[1 0 50]);
				allmetrics(di).fA_gain_tau(i) = ff.tau;
				cf = confint(ff);
				allmetrics(di).fA_gain_tau_err(i) = diff(cf(:,3));

				% fit the stimulus too
				S = MSGdata.PID(5e3+1:55e3,i);
				ff = fit(x(:),S,ft_stim,'StartPoint',[max(S) 1e2],'Lower',[0 50]);
				allmetrics(di).S_tau(i) = ff.tau;
			end
		end
	end
	save('.cache/tau_gain_MSG.mat','allmetrics')
end


figure('outerposition',[0 0 1250 602],'PaperUnits','points','PaperSize',[1250 602]); hold on
for di = 1

	% plot
	subplot(1,2,1); hold on
	for i = 1:max(allmetrics(di).orn)
		x = i;
		y = allmetrics(di).fA_tau(allmetrics(di).orn==i);
		y(y<51) = NaN; % ignore bad fits, because bound was hit
		y = 1e-3*nonnans(y);
		if length(y)>1
			plot(x,median(y),'k.','MarkerSize',40)
			plot(x,y,'k+')
		end
	end
	set(gca,'XLim',[.5 .5+max(allmetrics(di).orn)],'YLim',[0 50],'XTick',[1:max(allmetrics(di).orn)])
	xlabel('ORN')
	ylabel('\tau_{response} (s)')
	t = [orn_names{di} char(10) odour_names{di}];
	title(t)
	deintersectAxes;

	% plot
	subplot(1,2,2); hold on
	for i = 1:max(allmetrics(di).orn)
		x = i;
		y = allmetrics(di).fA_gain_tau(allmetrics(di).orn==i);
		y(y<51) = NaN; % ignore bad fits, because bound was hit
		y = nonnans(y);
		if length(y)>1
			plot(x,median(y),'k.','MarkerSize',40)
			plot(x,y,'k+')
		end
	end
	set(gca,'XLim',[.5 .5+max(allmetrics(di).orn)],'YLim',[0 1e3],'XTick',[1:max(allmetrics(di).orn)])
	xlabel('ORN')
	ylabel('\tau_{gain} (ms)')
	t = [orn_names{di} char(10) odour_names{di}];
	title(t)
	deintersectAxes;
end

prettyFig();

if being_published
	snapnow
	delete(gcf)
end



%% Version Info
%
pFooter;


