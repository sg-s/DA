
pHeader;


%% Estimating the timescale of gain control
% In this document, I attempt to determine the timescale of gain control from data where we present pulses of odorant on top of a background, at various time-points since background onset. 

% first, gather the data
root = '/data/DA-paper/data-for-paper/tau-gain/'; 
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

%
% look at kinetics

% figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on

% for si = 1:length(all_sensillum)
% 	ts = all_sensillum(si);

% 	subplot(1,length(all_sensillum),si); hold on
% 	for i = 1:length(all_pulse_times)
% 		tpt = all_pulse_times(i);
% 		a = (tpt+1)*1e3;
% 		z = a + .5e3;
% 		y = LFP(a:z,sensillum == ts & pulse_time == tpt);
% 		y = -mean(y,2); y = y - mean(y(1:10)); y = y/max(y); 
% 		plot(y,'Color',c(i,:))
% 	end
	

% end

% prettyFig();


%
% Now, I attempt to fit a NL model to this data. 

clear data
all_pulse_times = unique(pulse_time);
for i = 1:length(all_pulse_times)
	tpt = all_pulse_times(i);
	j = find(pulse_time == tpt & sensillum == 3,1,'first');
	data(i).stimulus = S(:,j);
	data(i).response = -LFP(:,j);
	%data(i).response(1:1e3) = NaN;
end

%% Version Info
%
pFooter;


