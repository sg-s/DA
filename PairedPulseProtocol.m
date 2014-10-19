% PairedPulseProtocol.m
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

% internal housekeeping: determine if being called by publish or not
calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end

%% Paired Pulse Protocol
% This document describes the paired pulse experiment. 

load('/local-data/DA-paper/ppp/2014_10_02_CSF2_EA_ab3_PairedPulses_2.mat')


%% Stimulus 
% The following figure shows the stimulus type we present: in one case, a small pulse is followed by a big pulse. In the other, a big pulse is followed by a small pulse. The panel on the top shows the two cases we consider. The panel on the bottom left shows one pair zoomed up, showing that the pulses are similar irrespective of order presented. In the final panel on the bottom right, we quantify the height of the pulses as a function of pulse number. We see that there is no significant difference between the blue and red curves, showing that the stimulus delivery is independent of the order of the pulses.  


figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
subplot(2,1,1), hold on
time = 1e-4:1e-4:1e-4*length(data(2).PID);
plot(time,mean2(data(2).PID))
plot(time,mean2(data(3).PID),'r')
xlabel('Time (s)')
ylabel('PID')

% zoom up on one
subplot(2,2,3), hold on
plot(time,mean2(data(2).PID))
plot(time,mean2(data(3).PID),'r')
set(gca,'XLim',[25.2 26.7])
xlabel('Time (s)')
ylabel('PID (V)')
legend({'Pp','pP'})


% variation of height with pulse #
subplot(2,2,4), hold on
npulses = sum(diff(ControlParadigm(2).Outputs(5,:)) == 1);
p_height = NaN(2,npulses);
P_height = NaN(2,npulses);
[bons,boffs]=ComputeOnsOffs(ControlParadigm(2).Outputs(6,:));
[fons,foffs]=ComputeOnsOffs(ControlParadigm(2).Outputs(5,:));
pulse_width = boffs(1)-bons(1);
for i = 1:npulses
	PID = mean2(data(2).PID);
	p_height(1,i)= max(PID(bons(i):boffs(i)+pulse_width/2));
	P_height(1,i)= max(PID(fons(i):foffs(i)+pulse_width/2));
end
[bons,boffs]=ComputeOnsOffs(ControlParadigm(3).Outputs(6,:));
[fons,foffs]=ComputeOnsOffs(ControlParadigm(3).Outputs(5,:));
pulse_width = boffs(1)-bons(1);
for i = 1:npulses
	PID = mean2(data(3).PID);
	p_height(2,i)= max(PID(bons(i):boffs(i)+pulse_width/2));
	P_height(2,i)= max(PID(fons(i):foffs(i)+pulse_width/2));
end
plot(p_height(1,:))
plot(p_height(2,:),'r')
plot(P_height(1,:))
plot(P_height(2,:),'r')
ylabel('Pulse Height (V)')
xlabel('Pulse #')


PrettyFig;

if being_published
	snapnow;
	delete(gcf)
end

%% ab3A Responses 
% In this section, we look at how the ab3A neuron responds to these pulses. 


% convert all the spiketimes in all the paradigms to a firing rate.
if ~exist('f')
	f(2).A = [];
	for i = 1:length(spikes)
		textbar(i,length(spikes))
		if ~isempty(spikes(i).A) && length(spikes(i).A) > 1
			[f_this,ft] = spiketimes2f(spikes(i).A);
			f(i).A  = mean2(f_this);
			f(i).time = ft;
		end
	end
end


figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
subplot(2,2,1:2), hold on
plot(f(2).time,f(2).A), hold on
plot(f(2).time,f(3).A,'r')
set(gca,'XLim',[25 38])
xlabel('Time (s)')
ylabel('Firing rate (Hz)')
legend({'Pp','pP'})
title('Lag=100ms')

% calcualte peak heights based on order 
Pp_paradigm = 2; pP_paradigm = 3;
small_peaks1 = zeros(1,10);
small_peaks2 = zeros(1,10);
large_peaks1 = zeros(1,10);
large_peaks2 = zeros(1,10);

% find when the small valve is on
small_pulse_1st = interp1(time,ControlParadigm(pP_paradigm).Outputs(6,:),f(pP_paradigm).time);
small_pulse_2nd = interp1(time,ControlParadigm(Pp_paradigm).Outputs(6,:),f(Pp_paradigm).time);
[ons_1st,offs_1st] = ComputeOnsOffs(small_pulse_1st);
[ons_2nd,offs_2nd] = ComputeOnsOffs(small_pulse_2nd);
dt=mean(diff(f(pP_paradigm).time));
% everything is delayed a 100ms, so correct for that
d =floor(.100/dt);
ons_1st = ons_1st + d; ons_2nd = ons_2nd + d; offs_1st = offs_1st + d; offs_2nd = offs_2nd + d;
for i = 1:length(small_peaks1)
	small_peaks1(i) = max(f(pP_paradigm).A(ons_1st(i):offs_1st(i)));
	small_peaks2(i) = max(f(Pp_paradigm).A(ons_2nd(i):offs_2nd(i)));
end

% find when the big valve is on
big_pulse_1st = interp1(time,ControlParadigm(Pp_paradigm).Outputs(5,:),f(Pp_paradigm).time);
big_pulse_2nd = interp1(time,ControlParadigm(pP_paradigm).Outputs(5,:),f(pP_paradigm).time);
[ons_1st,offs_1st] = ComputeOnsOffs(big_pulse_1st);
[ons_2nd,offs_2nd] = ComputeOnsOffs(big_pulse_2nd);
dt=mean(diff(f(pP_paradigm).time));
% everything is delayed a 100ms, so correct for that
d =floor(.100/dt);
ons_1st = ons_1st + d; ons_2nd = ons_2nd + d; offs_1st = offs_1st + d; offs_2nd = offs_2nd + d;
for i = 1:length(small_peaks1)
	large_peaks1(i) = max(f(Pp_paradigm).A(ons_1st(i):offs_1st(i)));
	large_peaks2(i) = max(f(pP_paradigm).A(ons_2nd(i):offs_2nd(i)));
end

subplot(2,2,3), hold on


plot(small_peaks1,small_peaks2,'g.','MarkerSize',35)
plot(large_peaks1,large_peaks2,'k.','MarkerSize',35)

legend({'small','Large'},'Location','northwest')
xlabel('Peak f when first (Hz)')
ylabel('Peak f when second (Hz)')
set(gca,'XLim',[20 120],'YLim',[20 120])
axis square

% plot a line of unity
plot([0 200],[0 200],'k--')

% plot lines to each of the clouds of points
fo=fit(small_peaks1(:),small_peaks2(:),'poly1');
plot(min(small_peaks1)-20:20+max(small_peaks1),fo(min(small_peaks1)-20:20+max(small_peaks1)),'g')


% plot lines to each of the clouds of points
fo=fit(large_peaks1(:),large_peaks2(:),'poly1');
plot(min(large_peaks1)-20:20+max(large_peaks1),fo(min(large_peaks1)-20:20+max(large_peaks1)),'k')

title('Lag=100ms')


% now calculate the ratios of each pair as a function of lag
r_small = NaN(10,5); r_large = NaN(10,5);
for i = 1:5

	Pp_paradigm = 2*i; % big pulse first paradigm
	pP_paradigm = 2*i + 1; % small pulse first paradigm
	time = 1e-4:1e-4:1e-4*length(data(Pp_paradigm).PID);

	% find when the big valve is on
	big_pulse_1st = interp1(time,ControlParadigm(Pp_paradigm).Outputs(5,:),f(Pp_paradigm).time);
	big_pulse_2nd = interp1(time,ControlParadigm(pP_paradigm).Outputs(5,:),f(pP_paradigm).time);
	[ons_1st,offs_1st] = ComputeOnsOffs(big_pulse_1st);
	[ons_2nd,offs_2nd] = ComputeOnsOffs(big_pulse_2nd);
	dt=mean(diff(f(pP_paradigm).time));
	% everything is delayed a 100ms, so correct for that
	d =floor(.100/dt);
	ons_1st = ons_1st + d; ons_2nd = ons_2nd + d; offs_1st = offs_1st + d; offs_2nd = offs_2nd + d;
	for j = 1:10
		large_peaks1(j) = max(f(Pp_paradigm).A(ons_1st(j):offs_1st(j)));
		large_peaks2(j) = max(f(pP_paradigm).A(ons_2nd(j):offs_2nd(j)));
	end

	r_large(:,i) = large_peaks1./large_peaks2;

	% find when the small valve is on
	small_pulse_1st = interp1(time,ControlParadigm(pP_paradigm).Outputs(6,:),f(pP_paradigm).time);
	small_pulse_2nd = interp1(time,ControlParadigm(Pp_paradigm).Outputs(6,:),f(Pp_paradigm).time);
	[ons_1st,offs_1st] = ComputeOnsOffs(small_pulse_1st);
	[ons_2nd,offs_2nd] = ComputeOnsOffs(small_pulse_2nd);
	dt=mean(diff(f(pP_paradigm).time));
	% everything is delayed a 100ms, so correct for that
	d =floor(.100/dt);
	ons_1st = ons_1st + d; ons_2nd = ons_2nd + d; offs_1st = offs_1st + d; offs_2nd = offs_2nd + d;
	for j = 1:10
		small_peaks1(j) = max(f(pP_paradigm).A(ons_1st(j):offs_1st(j)));
		small_peaks2(j) = max(f(Pp_paradigm).A(ons_2nd(j):offs_2nd(j)));
	end
	
	r_small(:,i) = small_peaks1./small_peaks2;

end


subplot(2,2,4), hold on
lag = 100:100:500;
plot(lag,r_small','g.','MarkerSize',30)
plot(lag,r_large','k.','MarkerSize',30)
xlabel('Lag between pulses (ms)')
ylabel('Ratio of pulses')

% plot a reference line
plot([0 max(lag)],[1 1],'k--')

PrettyFig;

if being_published
	snapnow;
	delete(gcf)
end



%% LN Model simulations 
% In this section, we use the LN model with parameters chosen from a representative dataset from Carlotta's flickering stimulus experiment to try to simulate what the neuron responses to these paired pulses will be. The filter is arbitrarily chosen from old data, and the output non-linearity is fit to the current data set. 

% load data
load('/local-data/DA-paper/data.mat')
td = 4;


% build a simple LN model
K = FindBestFilter(data(td).PID(500:end),data(td).ORN(500:end),[],'filter_length=201;');
filtertime = data(td).filtertime;



% now fit a nonlinearity to the data
K = K/1000;
load('/local-data/DA-paper/ppp/2014_10_02_CSF2_EA_ab3_PairedPulses_2.mat')
time = 1:length(data(2).PID); time = time*1e-4;
t = 0:3e-3:max(time);
PID = interp1(time,mean2(data(2).PID),t); PID(1) = 0;
fp2 = convolve(t,PID,K,filtertime);


xdata = fp2;
ydata = f(2).A;

% crop it to lose NaNs
ydata(isnan(xdata)) = [];
xdata(isnan(xdata)) = [];

xdata = xdata(:);
ydata = ydata(:);

fo=optimset('MaxFunEvals',1000,'Display','none');
x = lsqcurvefit(@hill,[max(ydata) 2 2],xdata,ydata,[max(ydata)/2 2 1],[2*max(ydata) max(ydata) 10],fo);
LNP2 = hill(x,fp2);

% now use this to predict the reversed order sequence
PID = interp1(time,mean2(data(3).PID),t); PID(1) = 0;
fp3 = convolve(t,PID,K,filtertime);
LNP3 = hill(x,fp3);

%%
% The following figure shows the output of the LN simulations using the same stimulus shown the ORN in the figure above. The LN model can capture the same qualitative behaviour we see in the neuron, indicating that the decrease in firing can be explained by a long inhibitory tail from a differentiating filter, and does not need modulation of gain. 

figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
subplot(2,2,1:2), hold on
plot(t,LNP2,'b'), hold on
plot(t,LNP3,'r'), hold on
set(gca,'XLim',[25 38])
xlabel('Time (s)')
ylabel('Firing rate (Hz)')
legend({'Pp','pP'})
title('LN Model Simulations Lag=100ms')


% calcualte peak heights based on order 
Pp_paradigm = 2; pP_paradigm = 3;
small_peaks1 = zeros(1,10);
small_peaks2 = zeros(1,10);
large_peaks1 = zeros(1,10);
large_peaks2 = zeros(1,10);

% find when the small valve is on
small_pulse_1st = interp1(time,ControlParadigm(pP_paradigm).Outputs(6,:),f(pP_paradigm).time);
small_pulse_2nd = interp1(time,ControlParadigm(Pp_paradigm).Outputs(6,:),f(Pp_paradigm).time);
[ons_1st,offs_1st] = ComputeOnsOffs(small_pulse_1st);
[ons_2nd,offs_2nd] = ComputeOnsOffs(small_pulse_2nd);
dt=3e-3;
% everything is delayed a 100ms, so correct for that
d =floor(.100/dt);
ons_1st = ons_1st + d; ons_2nd = ons_2nd + d; offs_1st = offs_1st + d; offs_2nd = offs_2nd + d;
for i = 1:length(small_peaks1)
	small_peaks1(i) = max(LNP3(ons_1st(i):offs_1st(i)));
	small_peaks2(i) = max(LNP2(ons_2nd(i):offs_2nd(i)));
end

% find when the big valve is on
big_pulse_1st = interp1(time,ControlParadigm(Pp_paradigm).Outputs(5,:),f(Pp_paradigm).time);
big_pulse_2nd = interp1(time,ControlParadigm(pP_paradigm).Outputs(5,:),f(pP_paradigm).time);
[ons_1st,offs_1st] = ComputeOnsOffs(big_pulse_1st);
[ons_2nd,offs_2nd] = ComputeOnsOffs(big_pulse_2nd);
dt=mean(diff(f(pP_paradigm).time));
% everything is delayed a 100ms, so correct for that
d =floor(.100/dt);
ons_1st = ons_1st + d; ons_2nd = ons_2nd + d; offs_1st = offs_1st + d; offs_2nd = offs_2nd + d;
for i = 1:length(small_peaks1)
	large_peaks1(i) = max(LNP2(ons_1st(i):offs_1st(i)));
	large_peaks2(i) = max(LNP3(ons_2nd(i):offs_2nd(i)));
end

subplot(2,2,3), hold on


plot(small_peaks1,small_peaks2,'g.','MarkerSize',35), hold on
plot(large_peaks1,large_peaks2,'k.','MarkerSize',35)

legend({'small','Large'},'Location','northwest')
xlabel('Peak f when first (Hz)')
ylabel('Peak f when second (Hz)')
set(gca,'XLim',[20 120],'YLim',[20 120])
axis square

% plot a line of unity
plot([0 200],[0 200],'k--')

% plot lines to each of the clouds of points
fo=fit(small_peaks1(:),small_peaks2(:),'poly1');
plot(min(small_peaks1)-20:20+max(small_peaks1),fo(min(small_peaks1)-20:20+max(small_peaks1)),'g')


% plot lines to each of the clouds of points
fo=fit(large_peaks1(:),large_peaks2(:),'poly1');
plot(min(large_peaks1)-20:20+max(large_peaks1),fo(min(large_peaks1)-20:20+max(large_peaks1)),'k')

title('Lag=100ms')


% now compute the LN model prediction for all lags.
fp = f;
for i = 2:11
	time = 1:length(data(i).PID); time = time*1e-4;
	t = 0:3e-3:max(time);
	PID = interp1(time,mean2(data(i).PID),t); PID(1) = 0;
	fp(i).A = convolve(t,PID,K,filtertime);
	fp(i).A = hill(x,fp(i).A);
end



% now calculate the ratios of each pair as a function of lag
r_small = NaN(10,5); r_large = NaN(10,5);
for i = 1:5

	Pp_paradigm = 2*i; % big pulse first paradigm
	pP_paradigm = 2*i + 1; % small pulse first paradigm
	time = 1e-4:1e-4:1e-4*length(data(Pp_paradigm).PID);

	% find when the big valve is on
	big_pulse_1st = interp1(time,ControlParadigm(Pp_paradigm).Outputs(5,:),f(Pp_paradigm).time);
	big_pulse_2nd = interp1(time,ControlParadigm(pP_paradigm).Outputs(5,:),f(pP_paradigm).time);
	[ons_1st,offs_1st] = ComputeOnsOffs(big_pulse_1st);
	[ons_2nd,offs_2nd] = ComputeOnsOffs(big_pulse_2nd);
	dt=mean(diff(f(pP_paradigm).time));
	% everything is delayed a 100ms, so correct for that
	d =floor(.100/dt);
	ons_1st = ons_1st + d; ons_2nd = ons_2nd + d; offs_1st = offs_1st + d; offs_2nd = offs_2nd + d;
	for j = 1:10
		large_peaks1(j) = max(fp(Pp_paradigm).A(ons_1st(j):offs_1st(j)));
		large_peaks2(j) = max(fp(pP_paradigm).A(ons_2nd(j):offs_2nd(j)));
	end

	r_large(:,i) = large_peaks1./large_peaks2;

	% find when the small valve is on
	small_pulse_1st = interp1(time,ControlParadigm(pP_paradigm).Outputs(6,:),f(pP_paradigm).time);
	small_pulse_2nd = interp1(time,ControlParadigm(Pp_paradigm).Outputs(6,:),f(Pp_paradigm).time);
	[ons_1st,offs_1st] = ComputeOnsOffs(small_pulse_1st);
	[ons_2nd,offs_2nd] = ComputeOnsOffs(small_pulse_2nd);
	dt=mean(diff(f(pP_paradigm).time));
	% everything is delayed a 100ms, so correct for that
	d =floor(.100/dt);
	ons_1st = ons_1st + d; ons_2nd = ons_2nd + d; offs_1st = offs_1st + d; offs_2nd = offs_2nd + d;
	for j = 1:10
		small_peaks1(j) = max(fp(pP_paradigm).A(ons_1st(j):offs_1st(j)));
		small_peaks2(j) = max(fp(Pp_paradigm).A(ons_2nd(j):offs_2nd(j)));
	end
	
	r_small(:,i) = small_peaks1./small_peaks2;

end


subplot(2,2,4), hold on
lag = 100:100:500;
plot(lag,r_small','g.','MarkerSize',30)
plot(lag,r_large','k.','MarkerSize',30)
xlabel('Lag between pulses (ms)')
ylabel('Ratio of pulses')

% plot a reference line
plot([0 max(lag)],[1 1],'k--')

PrettyFig;

if being_published
	snapnow;
	delete(gcf)
end

%% Gain Changes vs. Filter effects
% To clearly see if we can see cases where response decreases after a conditioning pulse due to the shape of the filter vs. response decreases due to a gain change, we present various pulses with and without conditioning pulses of different heights. 

load('/local-data/DA-paper/ppp2/ppp_data.mat')

% we organise the data in a 3D matrix. the first dimension indexes the trial # (globally), the second indexes what the conditioning pulse we use is (so 1 is no conditioning pulse), and the third dimension indexes the width of the pulse

% figure out what's in the data
conditioning_pulses  =unique([ppp_data.c_height]);
pulse_widths  =unique([ppp_data.width]);
probe_heights = unique([ppp_data.p_height]);

stim_height = NaN(300,length(conditioning_pulses),length(pulse_widths));
cond_stim_height = NaN(300,length(conditioning_pulses),length(pulse_widths));
resp_height = NaN(300,length(conditioning_pulses),length(pulse_widths));
cond_resp_height = NaN(300,length(conditioning_pulses),length(pulse_widths));
data_source = NaN(300,length(conditioning_pulses),length(pulse_widths));

% make matrices to store the time series averaged over all the pulses
dt = mean(diff(ppp_data(1).time));
before = floor(1/dt);
after = floor(2.5/dt);
triggered_data(length(probe_heights),length(conditioning_pulses),length(pulse_widths)).PID = zeros(before+after+1,1);
triggered_data(length(probe_heights),length(conditioning_pulses),length(pulse_widths)).ORN = zeros(before+after+1,1);

% c keeps a count, allowing us to fill the matrix correctly
c=ones(length(conditioning_pulses),length(pulse_widths));



for i = 1:length(ppp_data)
	% figure out what the conditioning pulse is 
	cp = find(ppp_data(i).c_height == conditioning_pulses);

	% figure out what the pulse width is
	pw = find(ppp_data(i).width == pulse_widths);

	% figure out what the probe pulse is 
	pp = find(ppp_data(i).p_height == probe_heights);

	if strcmp(ppp_data(i).neuron,'ab3') 
		% this is the neuron we want

		% estimate the baseline PID
		baseline = mean(ppp_data(i).PID(1:find(ppp_data(i).p_valve,1,'first')));
		thisPID = ppp_data(i).PID - baseline;

		% find valve ons and offs
		[ons,offs] = ComputeOnsOffs(ppp_data(i).p_valve);
		[ons2,offs2] = ComputeOnsOffs(ppp_data(i).c_valve);

		% look .1 second beyond the valve off 
		offs = offs + round(.1/dt);
		offs2 = offs2 + round(.1/dt);

		% % throw away 1st 5 pulses
		% ons(1:5) = [];
		% offs(1:5) = [];

		% if ~isempty(ons2)
		% 	ons2(1:5) = [];
		% 	offs2(1:5) = [];
		% end


		% extract values for each pulse
		for j = 1:length(ons)
			f_max = max(ppp_data(i).f(ons(j):offs(j)));
			p_max = max(thisPID(ons(j):offs(j)));

			if ~isempty(ons2)
				% also for the cond. pulse, if it exists
				f_max2 = max(ppp_data(i).f(ons2(j):offs2(j)));
				p_max2 = max(thisPID(ons2(j):offs2(j)));
				cond_stim_height(c(cp,pw),cp,pw) = p_max2;
				cond_resp_height(c(cp,pw),cp,pw) = f_max2;
			else
				cond_stim_height(c(cp,pw),cp,pw) = NaN;
				cond_resp_height(c(cp,pw),cp,pw) = NaN;
			end


			stim_height(c(cp,pw),cp,pw) = p_max;
			resp_height(c(cp,pw),cp,pw) = f_max;


			data_source(c(cp,pw),cp,pw) =  i;

			% add the full traces
			if isempty(triggered_data(pp,cp,pw).PID)
				triggered_data(pp,cp,pw).PID = ppp_data(i).PID(ons(j)-before:ons(j)+after)';
				triggered_data(pp,cp,pw).ORN = ppp_data(i).f(ons(j)-before:ons(j)+after);
			else
				triggered_data(pp,cp,pw).PID= [triggered_data(pp,cp,pw).PID ppp_data(i).PID(ons(j)-before:ons(j)+after)'];
				triggered_data(pp,cp,pw).ORN= [triggered_data(pp,cp,pw).ORN ppp_data(i).f(ons(j)-before:ons(j)+after)];
			end

			c(cp,pw) = c(cp,pw)+1;

		end
		
	end
end




%%
% The following figure shows response of the ab3 neuron to probe pulses of ethyl acetate, in the absence of conditioning, for pulse widths of 50ms. The panels on the left show the stimulus and the response for various amplitudes of the probe pulse, and the panel on the right shows the relationship between the peak of the stimulus and the peak of the response across all the data. Each dot on the right is a single pulse. Traces on the left are averaged over 5 trials. 



% make the plot
figure('outerposition',[0 0 1400 600],'PaperUnits','points','PaperSize',[1400 600]); hold on
a(1)=subplot(2,2,1); hold on

for i = 1:length(probe_heights)
	plot(a(1),dt*(-before:after),mean(triggered_data(i,1,1).PID'))
end
ylabel('Stimulus (a.u.)')

a(2) = subplot(2,2,3); hold on
for i = 1:length(probe_heights)
	plot(a(2),dt*(-before:after),mean(triggered_data(i,1,1).ORN'))
end
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')


a(3) = subplot(1,2,2); hold on
scatter(stim_height(:,1,1),resp_height(:,1,1))
ylabel('Peak Firing Rate (Hz)')
xlabel('Peak Stimulus (a.u.)')

PrettyFig;

if being_published
	snapnow;
	delete(gcf)
end


%%
% In comparison, the following figure shows the responses to a probe pulse following a conditioning pulse. The black dots are from the unconditioned responses, and the coloured dots are conditioning by a pulse just before the probe pulse (shown on the panels on the left)

this_cond = 2;

% make the plot
figure('outerposition',[0 0 1400 600],'PaperUnits','points','PaperSize',[1400 600]); hold on
a(1)=subplot(2,2,1); hold on


for i = 1:length(probe_heights)
	plot(a(1),dt*(-before:after),mean(triggered_data(i,this_cond,1).PID'))
end
ylabel('Stimulus (a.u.)')

a(2) = subplot(2,2,3); hold on

for i = 1:length(probe_heights)
	plot(a(2),dt*(-before:after),mean(triggered_data(i,this_cond,1).ORN'))
end
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')


a(3) = subplot(1,2,2); hold on
scatter(stim_height(:,1,1),resp_height(:,1,1),'k')
scatter(stim_height(:,this_cond,1),resp_height(:,this_cond,1))
ylabel('Peak Firing Rate (Hz)')
xlabel('Peak Stimulus (a.u.)')
set(a(3),'XLim',[0 3])

PrettyFig;

if being_published
	snapnow;
	delete(gcf)
end

%%
% What is the effect of the conditioning pulses on response properties of the neuron to the probe pulses? In the following figure, we plot the response vs. the stimulus for the different conditioning pulses used. Each dataset conditioned with a different pulse amplitude is shown in a different colour. 

% ignore large pulses
resp_height(stim_height(:,1,1)>2,1,1) = NaN;
stim_height(stim_height(:,1,1)>2,1,1) = NaN;

% make the plot
figure('outerposition',[0 0 1400 600],'PaperUnits','points','PaperSize',[1400 600]); hold on
a(1)=subplot(1,3,1); hold on

slopes = NaN*conditioning_pulses;
intercepts = NaN*conditioning_pulses;

col = jet(length(conditioning_pulses)-1);

for this_cond = 1:length(conditioning_pulses)-1
	x = stim_height(~isnan(stim_height(:,this_cond,1)),this_cond,1);
	y = resp_height(~isnan(resp_height(:,this_cond,1)),this_cond,1);
	scatter(x,y,32,col(this_cond,:))

	% fit a line
	ff = fit(x,y,'poly1');
	slopes(this_cond) = ff.p1;
	intercepts(this_cond) = ff.p2;

	% draw it
	plot(x,ff(x),'Color',col(this_cond,:))

end

a(2)=subplot(1,3,2); hold on
gain = slopes/slopes(1);

x = NaN*conditioning_pulses; x(1) = 0;
for i = 2:length(conditioning_pulses)-1
	x(i) =  mean2(cond_stim_height(:,i,1));
end
scatter(x,gain)
xlabel('Conditioning Pulse Height (a.u.)')
ylabel('Relative Gain to probe pulse')

a(3)=subplot(1,3,3); hold on
scatter(x,intercepts)
xlabel('Conditioning Pulse Height (a.u.)')
ylabel('Intercepts')
PrettyFig;

if being_published
	snapnow;
	delete(gcf)
end



% figure('outerposition',[0 0 1200 600],'PaperUnits','points','PaperSize',[1200 600]); hold on
% subplot(1,2,1), hold on

% for i = 1:length(conditioning_pulses)
% 	y = cond_stim_height(~isnan(cond_stim_height(:,i,1)),i,1);
% 	x = conditioning_pulses(i)*ones(1,length(y));
% 	scatter(x,y)
% end
% xlabel('Flow of conditioning pulse (mL/min)')
% ylabel('Stimulus peak (a.u.)')

% subplot(1,2,2), hold on

% for i = 1:length(conditioning_pulses)
% 	y = cond_resp_height(~isnan(cond_stim_height(:,i,1)),i,1);
% 	x = conditioning_pulses(i)*ones(1,length(y));
% 	scatter(x,y)
% end
% xlabel('Flow of conditioning pulse (mL/min)')
% ylabel('Response peak (Hz)')


%%
% How is the reponse to the probe stimulus changed by the conditioning pulse? 

