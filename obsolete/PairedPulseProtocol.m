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

load('/local-data/DA-paper/ppp1/2014_10_02_CSF2_EA_ab3_PairedPulses_2.mat')


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
load('/local-data/DA-paper/flickering-stim/data.mat')
td = 4;


% build a simple LN model
K = FindBestFilter(data(td).PID(500:end),data(td).ORN(500:end),[],'filter_length=201;');
filtertime = data(td).filtertime;



% now fit a nonlinearity to the data
K = K/1000;
load('/local-data/DA-paper/ppp1/2014_10_02_CSF2_EA_ab3_PairedPulses_2.mat')
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


% we're going to store all the metrics in long vectors, which each vector referring to one parameter
probe_type =  NaN(9000,1); % what was the intended probe?
stim_probe = NaN(9000,1); % stores the actual value of the probe
stim_cond = NaN(9000,1); % what was the actual coniditoning pulse?
cond_type = NaN(9000,1); % what was the intended conditioning stimulus?

width_type = NaN(9000,1); % what is the width? refer to pulse_widths to figure it out

resp_probe = NaN(9000,1); % response to probe
resp_cond = NaN(9000,1); % response to cond

data_source = NaN(9000,1); % where did this come from? 
neuron_id = NaN(9000,1);
lag = NaN(9000,1); % lag between pulses 


% make matrices to store the time series averaged over all the pulses
dt = mean(diff(ppp_data(1).time));
before = floor(1/dt);
after = floor(2.5/dt);
triggered_data(length(probe_heights),length(conditioning_pulses),length(pulse_widths)).PID = zeros(before+after+1,1);
triggered_data(length(probe_heights),length(conditioning_pulses),length(pulse_widths)).ORN = zeros(before+after+1,1);

% c keeps a count, allowing us to fill the matrix correctly
c=1;

allfiles = dir('/local-data/DA-paper/ppp2/raw/2014_10_09_CSF*.mat');

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

		% ignore first n pulse
		n=4;
		ons(1:n) = []; offs(1:n) =[]; 
		if ~isempty(ons2)
			ons2(1:n)= []; offs2(1:n) = [];
		end

		% extract values for each pulse
		for j = 1:length(ons)
			

			if ~isempty(ons2)
				
				%disp('probe and conditioning pulse')
			

				if ppp_data(i).width == 300
					% fall back to a simple algo
					[resp_cond(c),a] = max(ppp_data(i).f(ons2(j):offs2(j)));
					[resp_probe(c),b] = max(ppp_data(i).f(ons(j):offs(j)));
					[stim_cond(c),cc] = max(thisPID(ons2(j):offs2(j)));
					[stim_probe(c),d] = max(thisPID(ons(j)+15:offs(j)+15));

				else
					% find the peaks in the PID
					mpd = max([.1/dt (ppp_data(i).width/1000)/dt]);
					[~,loc]=findpeaks(thisPID(ons2(j):ons2(j)+floor(2/dt)),'SortStr','descend','MinPeakDistance',mpd);
					loc=loc(1:2); loc=loc+ons2(j)-1; loc=sort(loc,'ascend');
					stim_probe(c) = thisPID(loc(2));
					stim_cond(c) = thisPID(loc(1));		
					resp_cond(c) = max(ppp_data(i).f(loc(1):loc(2))); % maximum b/w PID max
					resp_probe(c) = max(ppp_data(i).f(loc(2):loc(2)+floor(1/dt)));	
				end
				%debug

				% figure, hold on
				% plot(thisPID)
				% scatter(loc(1),thisPID(loc(1)))
				% scatter(loc(2),thisPID(loc(2)))
				% title(i)
				% set(gca,'XLim',[ons2(j)-30 offs2(j)+200])

				% pause(1)

				% delete(gcf)

				

 				% pidlag = 20; 
				% ornlag = 80; % totally eyeballed. super inaccurate
				% f_max = max(ppp_data(i).f(ons(j)+ornlag:offs(j)+ornlag));
				% p_max = max(thisPID(ons(j)+pidlag:offs(j)+pidlag));
				% f_max2 = max(ppp_data(i).f(ons2(j)+ornlag:offs2(j)+ornlag));
				% p_max2 = max(thisPID(ons2(j)+pidlag:offs2(j)+pidlag));
				
				% resp_cond(c) = f_max2;
				% resp_probe(c) = f_max;
			else
				% disp('there is no conditioning pulse. just the probe')
				% so we compute the maximum from this on to the next (very generous, since we don't know the exact lag)
				a = ons(j);
				if length(ons) == j
					z = length(ppp_data(i).f);
				else
					z = ons(j+1);
				end
				f_max = max(ppp_data(i).f(a:z));
				p_max = max(thisPID(a:z));
				stim_probe(c) = p_max;
				stim_cond(c) = thisPID(a-10); % because there is no conditioning here, pick some point in the past
				resp_cond(c) = ppp_data(i).f(a-10);
				resp_probe(c) = f_max;

			end

			neuron_id(c) = find(strcmp(ppp_data(i).original_name,{allfiles.name}));
			if neuron_id(c) > 0
				neuron_id(c) = 1;
				% this is because the last four files are from the same neuron, according to Mahmut
			end
			data_source(c) = i;



			% add the full traces
			if isempty(triggered_data(pp,cp,pw).PID)
				triggered_data(pp,cp,pw).PID = ppp_data(i).PID(ons(j)-before:ons(j)+after)';
				triggered_data(pp,cp,pw).ORN = ppp_data(i).f(ons(j)-before:ons(j)+after);
			else
				triggered_data(pp,cp,pw).PID= [triggered_data(pp,cp,pw).PID ppp_data(i).PID(ons(j)-before:ons(j)+after)'];
				triggered_data(pp,cp,pw).ORN= [triggered_data(pp,cp,pw).ORN ppp_data(i).f(ons(j)-before:ons(j)+after)];
			end


			width_type(c) = pw;
			probe_type(c) = pp;
			cond_type(c) = cp;
			c = c+1;

		end
		
	end
end




%%
% The following figure shows response of the ab3 neuron to probe pulses of ethyl acetate, in the absence of conditioning, for pulse widths of 50ms. The panels on the left show the stimulus and the response for various amplitudes of the probe pulse, and the panel on the right shows the relationship between the peak of the stimulus and the peak of the response across all the data. Each dot on the right is a single pulse. Traces on the left are averaged over 5 trials. 



% make the plot
figure('outerposition',[0 0 1400 600],'PaperUnits','points','PaperSize',[1400 600]); hold on
a(1)=subplot(2,2,1); hold on

col = pmkmp(length(probe_heights),'IsoL');
for i = 1:length(probe_heights)
	plot(a(1),dt*(-before:after),mean(triggered_data(i,1,1).PID'),'Color',col(i,:))
end
ylabel('Stimulus (a.u.)')

a(2) = subplot(2,2,3); hold on
for i = 1:length(probe_heights)
	plot(a(2),dt*(-before:after),mean(triggered_data(i,1,1).ORN'),'Color',col(i,:))
end
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')


a(3) = subplot(1,2,2); hold on
x = stim_probe(cond_type==1 & width_type == 1);
y = resp_probe(cond_type==1 & width_type == 1);
scatter(x,y,32)
ylabel('Peak Firing Rate (Hz)')
xlabel('Peak Stimulus (a.u.)')

PrettyFig;

if being_published
	snapnow;
	delete(gcf)
end


%%
% What is the effect of a conditioning pulse on the response of the neuron? In the following figure, we compare the response of the neuron to a probe pulse with and without a conditioning pulse just before the probe pulse. 

this_probe = 3;
this_width = 1;

% make the plot
figure('outerposition',[0 0 1400 600],'PaperUnits','points','PaperSize',[1400 600]); hold on
a(1)=subplot(2,2,1); hold on

col = pmkmp(length(conditioning_pulses),'IsoL');
for i = 1:length(conditioning_pulses)
	plot(a(1),dt*(-before:after),mean(triggered_data(this_probe,i,this_width).PID'),'Color',col(i,:))
end
ylabel('Stimulus (a.u.)')

a(2) = subplot(2,2,3); hold on

for i = 1:length(conditioning_pulses)
	plot(a(2),dt*(-before:after),mean(triggered_data(this_probe,i,this_width).ORN'),'Color',col(i,:))
end
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')


a(3) = subplot(2,2,2); hold on
x = conditioning_pulses(cond_type(probe_type==this_probe & width_type == this_width));
y = stim_cond(probe_type==this_probe & width_type == this_width);
scatter(x,y,32)
xlabel('Conditioning Pulse Flow rate (mL/min)')
ylabel('Cond. Pulse (a.u.)')

a(4) = subplot(2,2,4); hold on
x = stim_cond(probe_type==this_probe & width_type == this_width);
y = resp_probe(probe_type==this_probe & width_type == this_width);
scatter(x,y,32)

ylabel('Response to Probe (Hz)')
xlabel('Conditioning pulse peak (a.u.)')


PrettyFig;

if being_published
	snapnow;
	delete(gcf)
end


%% 
% What is the effect of the conditioning pulse on the response properties of the neuron? In the following analysis, we start by pulling out the responses of each neuron to all probe pulses in the absence of any conditioning for the 50ms probes.

% ignore super large probes
resp_probe(stim_probe>5) = NaN;
stim_probe(stim_probe>5) = NaN;

this_width = 1;
figure('outerposition',[0 0 1200 600],'PaperUnits','points','PaperSize',[1200 600]); hold on
subplot(1,2,1), hold on
neurons = unique(neuron_id(~isnan(neuron_id)));
n = length(neurons);
probe_mins = NaN(n,1);
probe_maxs = NaN(n,1);
for i  = 1:length(neurons)
	this_neuron = neurons(i);
	x = stim_probe(neuron_id==this_neuron & width_type == this_width & cond_type == 1);
	y = resp_probe(neuron_id==this_neuron & width_type == this_width & cond_type == 1);
	x = x(~isnan(x)); y= y(~isnan(y));
	scatter(x,y,'k')
	if ~isempty(x)
		ff = fit(x,y,'poly1');
		plot(x,ff(x),'k')
		probe_mins(i) = min(x);
		probe_maxs(i) = max(x);
	end
	
	
end
set(gca,'XLim',[0 3],'YLim',[0 160])
xlabel('Probe Stimulus Height (V)')
ylabel('Probe Response Height (Hz)')
title('Responses to probe alone grouped by neuron')

PrettyFig;

if being_published
	snapnow;
end


%% 
% We then show the responses of this neuron to a probe following a conditioning pulse. 

probe_ranges = probe_maxs - probe_mins;
tolerance = .1; % tolerance for probe ranges
col = pmkmp(length(conditioning_pulses)-1,'CubicL');
col = [0 0 0; col];
for i  = 1:length(neurons)
	for j = 2:length(conditioning_pulses)
		this_cond_probe_max = max((stim_probe(neuron_id==i & width_type == this_width & cond_type == j)));
		this_cond_probe_min = min((stim_probe(neuron_id==i & width_type == this_width & cond_type == j)));
		if ~isempty(this_cond_probe_min)
			% check if the range of the probes are OK

			%if abs((this_cond_probe_max - this_cond_probe_min) - probe_ranges(i)) < tolerance*probe_ranges(i)
				% yay! this is OK
			
				x = stim_probe(neuron_id==i & width_type == this_width & cond_type == j);
				y = resp_probe(neuron_id==i & width_type == this_width & cond_type == j);
				scatter(x,y,64,col(j,:))
				x = x(~isnan(x)); y= y(~isnan(y));
				if ~isempty(x)
					ff = fit(x,y,'poly1');
					plot(x,ff(x),'Color',col(j,:))

				end
				
			% else
			% 	(this_cond_probe_max - this_cond_probe_min) 
			% end
		end

	end
	clear j
end
clear i
title('Unconditioned (black) vs conditioned responses (coloured)')

% make labels
cond_value = [];
for i = 1:length(conditioning_pulses)
	cond_value(i) = mean(stim_cond(width_type==this_width & cond_type==i));
end
cond_value = reshape(repmat(cond_value,2,1),1,length(conditioning_pulses)*2);
for i = 1:length(cond_value)
	cond_label{i} = (oval(cond_value(i),2));
end
legend(cond_label,'Location','southeast')


PrettyFig;

if being_published
	snapnow;
end


%%
% Now, we perform a similar analysis on the 300ms pulses (panel on the right)

title('50ms')
this_width = 2;

subplot(1,2,2), hold on
neurons = unique(neuron_id(~isnan(neuron_id)));
n = length(neurons);
probe_mins = NaN(n,1);
probe_maxs = NaN(n,1);
for i  = 1:length(neurons)
	this_neuron = neurons(i);
	x = stim_probe(neuron_id==this_neuron & width_type == this_width & cond_type == 1);
	y = resp_probe(neuron_id==this_neuron & width_type == this_width & cond_type == 1);
	x = x(~isnan(x)); y= y(~isnan(y));
	scatter(x,y,'k')
	if ~isempty(x)
		ff = fit(x,y,'poly1');
		plot(x,ff(x),'k')
		probe_mins(i) = min(x);
		probe_maxs(i) = max(x);
	end
	
	
end
set(gca,'XLim',[0 3],'YLim',[0 160])
xlabel('Probe Stimulus Height (V)')
ylabel('Probe Response Height (Hz)')



probe_ranges = probe_maxs - probe_mins;
tolerance = .1; % tolerance for probe ranges
col = pmkmp(length(conditioning_pulses)-1,'CubicL');
col = [0 0 0; col];
for i  = 1:length(neurons)
	for j = 2:length(conditioning_pulses)
		this_cond_probe_max = max((stim_probe(neuron_id==i & width_type == this_width & cond_type == j)));
		this_cond_probe_min = min((stim_probe(neuron_id==i & width_type == this_width & cond_type == j)));
		if ~isempty(this_cond_probe_min)
			% check if the range of the probes are OK

			%if abs((this_cond_probe_max - this_cond_probe_min) - probe_ranges(i)) < tolerance*probe_ranges(i)
				% yay! this is OK
			
				x = stim_probe(neuron_id==i & width_type == this_width & cond_type == j);
				y = resp_probe(neuron_id==i & width_type == this_width & cond_type == j);
				scatter(x,y,64,col(j,:))
				x = x(~isnan(x)); y= y(~isnan(y));
				if ~isempty(x)
					ff = fit(x,y,'poly1');
					plot(x,ff(x),'Color',col(j,:))

				end
				
			% else
			% 	(this_cond_probe_max - this_cond_probe_min) 
			% end
		end

	end
	clear j
end
clear i
title('300ms pulses')

cond_value = [];
% make labels
for i = 1:length(conditioning_pulses)
	cond_value(i) = mean(stim_cond(width_type==this_width & cond_type==i));
end
cond_value = reshape(repmat(cond_value,2,1),1,length(conditioning_pulses)*2);
for i = 1:length(cond_value)
	cond_label{i} = (oval(cond_value(i),2));
end
legend(cond_label,'Location','southeast')

PrettyFig;


return

%%
% What is the effect of the conditioning pulses on response properties of the neuron to the probe pulses? In the following figure, we plot the response vs. the stimulus for the different conditioning pulses used. Each dataset conditioned with a different pulse amplitude is shown in a different colour. 

this_width = 1;

% ignore large pulses
resp_probe(stim_probe>2) = NaN;
stim_probe(stim_probe>2) = NaN;

% make the plot
figure('outerposition',[0 0 1200 600],'PaperUnits','points','PaperSize',[1200 600]); hold on
a(1)=subplot(1,2,1); hold on

slopes = NaN*conditioning_pulses;
intercepts = NaN*conditioning_pulses;

col = pmkmp(length(conditioning_pulses),'LinLhot');

for this_cond = 1:length(conditioning_pulses)
	x = stim_probe(cond_type==this_cond & width_type == this_width);
	y = resp_probe(cond_type==this_cond & width_type == this_width);
	scatter(x,y,32,col(this_cond,:))

	x = x(~isnan(x)); y = y(~isnan(y));

	% fit a line
	ff = fit(x,y,'poly1');
	slopes(this_cond) = ff.p1;
	intercepts(this_cond) = ff.p2;

	% draw it
	plot(x,ff(x),'Color',col(this_cond,:))

end
xlabel('Stimulus Probe Peak (a.u.)')
ylabel('Peak Response to Probe (Hz)')

a(2)=subplot(1,2,2); hold on
gain = slopes/slopes(1);

x = NaN*conditioning_pulses; x(1) = 0;
for i = 2:length(conditioning_pulses)
	x(i) =  mean(stim_cond(cond_type==i & width_type == this_width));
end

scatter(x,gain,164,col,'filled')
xlabel('Conditioning Pulse Height (a.u.)')
ylabel('Relative Gain to probe pulse')

PrettyFig;

if being_published
	snapnow;
	delete(gcf)
end

%%
% We now repeat the same analysis for the 300ms pulses. 


this_width = 2;

% ignore large pulses
resp_probe(stim_probe>2) = NaN;
stim_probe(stim_probe>2) = NaN;

% make the plot
figure('outerposition',[0 0 1200 600],'PaperUnits','points','PaperSize',[1200 600]); hold on
a(1)=subplot(1,2,1); hold on

slopes = NaN*conditioning_pulses;
intercepts = NaN*conditioning_pulses;

col = pmkmp(length(conditioning_pulses),'LinLhot');

for this_cond = 1:length(conditioning_pulses)
	x = stim_probe(cond_type==this_cond & width_type == this_width);
	y = resp_probe(cond_type==this_cond & width_type == this_width);
	scatter(x,y,32,col(this_cond,:))

	x = x(~isnan(x)); y = y(~isnan(y));

	if ~isempty(x)
		% fit a line
		ff = fit(x,y,'poly1');
		slopes(this_cond) = ff.p1;
		intercepts(this_cond) = ff.p2;

		% draw it
		plot(x,ff(x),'Color',col(this_cond,:))
	end

	

end
xlabel('Stimulus Probe Peak (a.u.)')
ylabel('Peak Response to Probe (Hz)')

a(2)=subplot(1,2,2); hold on
gain = slopes/slopes(1);

x = NaN*conditioning_pulses; x(1) = 0;
for i = 2:length(conditioning_pulses)
	x(i) =  mean(stim_cond(cond_type==i & width_type == this_width));
end

scatter(x,gain,164,col,'filled')
xlabel('Conditioning Pulse Height (a.u.)')
ylabel('Relative Gain to probe pulse')

PrettyFig;

if being_published
	snapnow;
	delete(gcf)
end
