
pHeader;

%% Odour plume statistics
% In this document, I look at the statistics of real odour plumes. I generate them using a fan blowing over a vial of either Apple Cider Vinegar or ethyl acetate. 

% get data 
root = getPath(dataManager,'f24ec4a8dd9169efbfab300baa193579');

% first, gather the data
allfiles = dir([root '/*.kontroller']);
allfiles(cellfun(@(x) strcmp(x(1),'.'), {allfiles.name})) = [];

odour_x_location = [];
odour_y_location = [];
fan_power = [];
PID_location = [];
odour_name = {};
PID = [];


for i = 1:length(allfiles)
	clear data spikes ControlParadigm

	load([root '/' allfiles(i).name],'-mat')
	%disp({ControlParadigm.Name})

	a = strfind(allfiles(i).name,'PID_')+4;
	z = strfind(allfiles(i).name,'cm')-1;
	this_PID_location = str2double(allfiles(i).name(a:z));

	this_PID = vertcat(data.PID);
	this_PID = this_PID';
	this_PID = this_PID(1:10:end,:);

	z = strfind(allfiles(i).name,'_fan')-1;
	this_odour = allfiles(i).name(12:z);

	this_odour = repmat({this_odour},size(this_PID,2),1);
	this_PID_location = repmat(this_PID_location,size(this_PID,2),1);

	% figure out the positions of the odour for each trial
	n_trials_per_paradigm = cellfun(@width,{data.PID});
	this_x = [];
	this_y = [];
	this_fan_power = [];
	for j = 1:length(data)
		%disp(ControlParadigm(j).Name)
		z = strfind(ControlParadigm(j).Name,'cm'); 
		if length(z) > 1
			yz = z(2);
			z = z(1);
			a = strfind(ControlParadigm(j).Name,'_'); ya = a(end); a = a(1);
			x = str2double(ControlParadigm(j).Name(a+1:z-1));
			y = str2double(ControlParadigm(j).Name(ya+1:yz-1));
			this_x = [this_x; zeros(n_trials_per_paradigm(j),1)+x];
			this_y = [this_y; zeros(n_trials_per_paradigm(j),1)+y];
		else
			a = strfind(ControlParadigm(j).Name,'_');
			x = str2double(ControlParadigm(j).Name(a+1:z-1));
			this_x = [this_x; zeros(n_trials_per_paradigm(j),1)+x];
			this_y = [this_y; zeros(n_trials_per_paradigm(j),1)];
		end

		% figure out fan power
		if isempty(strfind(ControlParadigm(j).Name,'fan_off'))
			this_fan_power = [this_fan_power; ones(n_trials_per_paradigm(j),1)];
		else
			this_fan_power = [this_fan_power; zeros(n_trials_per_paradigm(j),1)];
		end
		
	end


	% consolidate
	PID = [PID this_PID];
	PID_location = [PID_location; this_PID_location];
	odour_x_location = [odour_x_location; this_x];
	odour_y_location = [odour_y_location; this_y];
	fan_power = [fan_power; this_fan_power];
	odour_name = [odour_name; this_odour];

end

clearvars -except PID odour_name PID_location odour_x_location odour_y_location being_published fan_power

% this defines the noise floor in the PID measurements 
noise_floor = 0.0036;

%%
% In the following figure, I show some sample traces of PID measurements of Apple Cider Vinegar plumes, measured at 35 cm from the fan. I plot the data for various values of the distance of the odour source from the fan. Note that odour at 35 cm means that the odour is directly below the PID inlet tube. Red traces are obtained with the fan turned off. I also plot the stimulus distributions of these time series in the final plot. The dashed line in the third plot indicates the noise floor of the PID measurement. 

figure('outerposition',[0 0 1502 501],'PaperUnits','points','PaperSize',[1502 501]); hold on

all_x_locations = unique(odour_x_location(strcmp(odour_name,'ACV')));
all_y_locations = unique(odour_y_location(strcmp(odour_name,'ACV')));
time = 1e-3*(1:length(PID));

c = 1;
cc = parula(3);
for i = 1:length(all_x_locations)
	for j = 1:length(all_y_locations)
		plot_these = find(odour_x_location == all_x_locations(i) & odour_y_location == all_y_locations(j) & strcmp(odour_name,'ACV') & fan_power);
		subplot(1,3,c); hold on
		ylabel('Stimulus (V)')
		xlabel('Time (s)')
		plot(time,PID(:,plot_these(end)),'Color',cc(c,:))
		t = ['x=' oval(odour_x_location(plot_these(1))) ', y=' oval(odour_y_location(plot_these(1)))];
		title(t)

		% also plot the statistics
		y = PID(:,plot_these); y = y(:); sy = sort(y); y = y - mean(sy(1:1e3));
		[hy,hx] = histcounts(y,100); hy = hy/sum(hy);
		hx = hx(1:end-1) + mean(diff(hx));
		subplot(1,3,3); hold on
		plot(hx,hy,'+','Color',cc(c,:))
		set(gca,'XScale','log','YScale','log','XLim',[1e-4 10],'XTick',logspace(-4,1,6))
		xlabel('Stimulus (V)')
		ylabel('Probability')
		
		
		plot_these = find(odour_x_location == all_x_locations(i) & odour_y_location == all_y_locations(j) & strcmp(odour_name,'ACV') & ~fan_power);
		if ~isempty(plot_these)
			subplot(1,3,c); hold on
			plot(time,PID(:,plot_these(1)),'r')
			y = PID(:,plot_these); y = y(:); sy = sort(y); y = y - mean(sy(1:1e3));
			[hy,hx] = histcounts(y,100); hy = hy/sum(hy);
			hx = hx(1:end-1) + mean(diff(hx));
			subplot(1,3,3); hold on
			plot(hx,hy,'+','Color','r')
		end

		c = c+1;

	end
end

% plot the noise floor
subplot(1,3,3); hold on
plot([noise_floor noise_floor], [1e-6 1],'k--')

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%%
% The PID could barely detect apple cider vinegar. I therefore switched to ethyl acetate as a proxy and made the same plots as before.  

figure('outerposition',[0 0 1230 901],'PaperUnits','points','PaperSize',[1230 901]); hold on

all_x_locations = unique(odour_x_location(strcmp(odour_name,'ethyl_acetate')));
all_y_locations = unique(odour_y_location(strcmp(odour_name,'ethyl_acetate')));
time = 1e-3*(1:length(PID));


c = 1;
cc = parula(6);
for i = 1:length(all_x_locations)
	for j = 1:length(all_y_locations)
		plot_these = find(odour_x_location == all_x_locations(i) & odour_y_location == all_y_locations(j) & strcmp(odour_name,'ethyl_acetate') & fan_power);
		if isempty(plot_these)
		else
			subplot(2,3,c); hold on
			y = PID(:,plot_these(1)); sy = sort(y); y = y - mean(sy(1:1e3));
			plot(time,y,'Color',cc(c,:))
			set(gca,'XLim',[10 15],'YLim',[-0.1 3],'YScale','linear')
			t = ['x=' oval(odour_x_location(plot_these(1))) ', y=' oval(odour_y_location(plot_these(1)))];
			title(t)

			% also plot the statistics
			y = PID(:,plot_these); y = y(:); sy = sort(y); y = y - mean(sy(1:1e3));
			[hy,hx] = histcounts(y,100); hy = hy/sum(hy);
			hx = hx(1:end-1) + mean(diff(hx));
			subplot(2,3,6); hold on
			plot(hx,hy,'+','Color',cc(c,:))
			set(gca,'XScale','log','YScale','log','XLim',[1e-4 10],'XTick',logspace(-4,1,6))
			xlabel('Stimulus (V)')
			ylabel('Probability')
			c = c+1;
		end
	end
end

% plot the noise floor
subplot(2,3,6); hold on
plot([noise_floor noise_floor], [1e-6 1],'k--')


prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end


%% Neuron responses
% In this section I plot response of neurons to these "real" odor plumes. We placed the PID downstream of the neuron, and the odour was about 10cm from the fly (the closest we could get). The fan was around 30 cm away, and the speed of the fan was varied in these recordings (using varying duty cycles @100Hz to the power system). 


odor = {};
ab = [];
sensillum = [];
S = [];
fA = [];
LFP = [];
fly = [];
odor_position = {};


% assemble data
root = [getPath(dataManager,'4c7dcd527caa5d3a3f6601aa6139e639') oss];
allfiles = dir([root '*.kontroller']);
allfiles(cellfun(@(x) strcmp(x(1),'.'), {allfiles.name})) = [];
for i = 1:length(allfiles)
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


	% figure out the odour
	this_odor = '';
	if any(strfind(allfiles(i).name,'2-butanone'))
		this_odor = '2-butanone';
	elseif any(strfind(allfiles(i).name,'ACV'))
		this_odor = 'ACV';
	elseif any(strfind(allfiles(i).name,'2ac'))
		this_odor = '2ac';
	elseif any(strfind(allfiles(i).name,'ethyl-butyrate'))
		this_odor = 'ethyl-butyrate';
	else 
		error('cant match odor')
	end
	this_odor = repmat({this_odor},length(this_ab),1);

	% figure out the odor position
	this_op = 'center';
	if any(strfind(allfiles(i).name,'left'))
		this_op = 'left';
	end
	this_op = repmat({this_op},length(this_ab),1);

	% consolidate
	fly = [fly; this_fly];
	odor = [odor; this_odor];
	odor_position = [odor_position; this_op];
	sensillum = [sensillum; this_sensillum];
	S = [S this_S];
	fA = [fA this_fA];
	LFP = [LFP this_LFP];
	ab = [ab; this_ab];

end

% remove bad LFP trials
rm_this = isnan(sum(LFP));
odor(rm_this) = [];
ab(rm_this) = [];
sensillum(rm_this) = [];
S(:,rm_this) = [];
fA(:,rm_this) = [];
LFP(:,rm_this) = [];
fly(rm_this) = [];
odor_position(rm_this) = [];

% remove baseline from LFP, and minimum from S
S = S - min(min(S));
for i = 1:size(LFP,2)
	LFP(:,i) = LFP(:,i) - LFP(1,i);
end

LFP = LFP*10;

all_odors = unique(odor);
c = 0;
ncols = 7;
time = 1e-3*(1:length(S));

%%
% In the following figure, I plot the stimulus, LFP and the firing rate for all the neurons we recorded from to get a sense of the data. 

figure('outerposition',[0 0 1411 900],'PaperUnits','points','PaperSize',[1411 900]); hold on
for i = 1:length(all_odors)
	this_odor = all_odors{i};
	for this_ab = 2:3
		if any(find(strcmp(odor,this_odor) & ab == this_ab))
			c = c+1;
			
			these = find(strcmp(odor,this_odor) & ab == this_ab);
			cc = parula(length(these)+1);
			for j = 1:length(these)
				this = these(j);
				subplot(3,ncols,c); hold on
				plot(time,S(:,this),'Color',cc(j,:))
				title(['ab' oval(this_ab) ' ' this_odor])

				subplot(3,ncols,ncols+c); hold on
				plot(time,LFP(:,this),'Color',cc(j,:))

				subplot(3,ncols,ncols*2+c); hold on
				plot(time,fA(:,this),'Color',cc(j,:))
			end
		end
	end
end


for i = 1:ncols
	subplot(3,ncols,i);
	set(gca,'XLim',[15 20],'YLim',[0 1])

	subplot(3,ncols,ncols+i);
	set(gca,'XLim',[15 20],'YLim',[-12 2])

	subplot(3,ncols,2*ncols+i);
	set(gca,'XLim',[15 20],'YLim',[0 120])
end

prettyFig();

if being_published
	snapnow
	delete(gcf)
end


%%
% Now I plot the distributions of the stimulus, the LFP and the response for every trial. 


all_odors = unique(odor);
c = 0;
ncols = 7;
time = 1e-3*(1:length(S));

figure('outerposition',[0 0 1411 900],'PaperUnits','points','PaperSize',[1411 900]); hold on
for i = 1:length(all_odors)
	this_odor = all_odors{i};
	for this_ab = 2:3
		if any(find(strcmp(odor,this_odor) & ab == this_ab))
			c = c+1;
			these = find(strcmp(odor,this_odor) & ab == this_ab);
			cc = parula(length(these)+1);
			for j = 1:length(these)
				this = these(j);
				subplot(3,ncols,c); hold on
				[hy,hx] = histcounts(S(:,this),min([100 length(unique(S(:,this)))]));
				hx = hx(1:end-1) + mean(diff(hx)); hy = hy/sum(hy);
				plot(hx,hy,'Color',cc(j,:))
				title(['ab' oval(this_ab) ' ' this_odor])

				subplot(3,ncols,ncols+c); hold on
				[hy,hx] = histcounts(LFP(:,this),100);
				hx = hx(1:end-1) + mean(diff(hx)); hy = hy/sum(hy);
				plot(hx,hy,'Color',cc(j,:))

				subplot(3,ncols,2*ncols+c); hold on
				if max(fA(:,this)) > 0
					[hy,hx] = histcounts(fA(:,this),100);
					hx = hx(1:end-1) + mean(diff(hx)); hy = hy/sum(hy);
					plot(hx,hy,'Color',cc(j,:))
				end
			end
		end
	end
end

for i = 1:ncols
	subplot(3,ncols,i);
	set(gca,'XScale','log')

	% subplot(3,ncols,ncols+i);
	% set(gca,'XLim',[15 25],'YLim',[-12 2])

	subplot(3,ncols,2*ncols+i);
	set(gca,'XScale','log','XLim',[1 200],'XTick',[1 10 100])
end

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% Now, I plot the cross-correlation functions between the stimulus and the LFP and the firing rate for all the data to get a feel for how well the measured stimulus can predict the observed response. 

all_odors = unique(odor);
c = 0;
ncols = 7;
time = 1e-3*(1:length(S));
peak_xcorr_LFP = NaN*ab;
peak_xcorr_fA = NaN*ab;

figure('outerposition',[0 0 1411 600],'PaperUnits','points','PaperSize',[1411 600]); hold on
for i = 1:length(all_odors)
	this_odor = all_odors{i};
	for this_ab = 2:3
		if any(find(strcmp(odor,this_odor) & ab == this_ab))
			c = c+1;
			these = find(strcmp(odor,this_odor) & ab == this_ab);
			cc = parula(length(these)+1);
			for j = 1:length(these)
				this = these(j);
				

				s = reshape(S(:,this),1e3,60); 
				x = reshape(LFP(:,this),1e3,60); 
				f = reshape(fA(:,this),1e3,60); 

				xcorr_x = NaN(2e3-1,60);
				xcorr_f = NaN(2e3-1,60);

				for k = 1:60
					s(:,k) = s(:,k) - mean(s(:,k));
					f(:,k) = f(:,k) - mean(f(:,k));
					x(:,k) = x(:,k) - mean(x(:,k));

					s(:,k) = s(:,k)/std(s(:,k));
					f(:,k) = f(:,k)/std(f(:,k));
					x(:,k) = x(:,k)/std(x(:,k));

					xcorr_x(:,k) = xcorr(x(:,k),s(:,k));
					if max(f(:,k)) > 0
						xcorr_f(:,k) = xcorr(f(:,k),s(:,k));
					end
				end
				xcorr_f = nanmean(xcorr_f,2);
				xcorr_x = nanmean(xcorr_x,2);

				xcorr_x = xcorr_x/1e3;
				xcorr_f = xcorr_f/1e3;

				peak_xcorr_LFP(this) = max(abs(xcorr_x));
				peak_xcorr_fA(this) = max(abs(xcorr_f));

				subplot(2,ncols,c); hold on
				lags = 1e-3*(1:length(xcorr_f)) - 1;
				plot(lags,xcorr_x,'Color',cc(j,:));
				title(['ab' oval(this_ab) ' ' this_odor])

				subplot(2,ncols,ncols+c); hold on
				plot(lags,xcorr_f,'Color',cc(j,:));
				
			end
		end
	end
end

for i = 1:ncols
	subplot(2,ncols,i);
	set(gca,'YLim',[-1 .5])
	xlabel('Lag (s)')
	ylabel('S \otimes LFP')


	subplot(2,ncols,ncols+i);
	set(gca,'YLim',[-.5 1])
	xlabel('Lag (s)')
	ylabel('S \otimes F')
end

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% We see that the overall cross correlation is pretty poor. The best case seems to be with 2-butanone (that the PID is extremely sensitive to). I now back out filters and plot the projected stimulus vs. the response for these cases. 

do_these = find(peak_xcorr_fA>.4);

figure('outerposition',[0 0 1501 703],'PaperUnits','points','PaperSize',[1501 703]); hold on
for i = 1:length(do_these)
	do_this = do_these(i);
	K = fitFilter2Data(S(:,do_this),fA(:,do_this),'offset',200,'reg',1);
	K = K(100:end-100);
	filtertime = 1e-3*(1:length(K)) - .1;
	subplot(2,length(do_these),i); hold on
	plot(filtertime,K*1e3,'k')
	title(['ab' oval(ab(do_this)) ' ' odor{do_this}])
	ylabel('Filter (a.u.)')
	xlabel('Lag (s)')

	fp = convolve(time,S(:,do_this),K,filtertime);
	subplot(2,length(do_these),length(do_these)+i); hold on
	plot(fp,fA(:,do_this),'k')
	xlabel('Proj. Stim. (V)')
	ylabel('Firing rate (Hz)')

end

prettyFig();

if being_published
	snapnow
	delete(gcf)
end


%% Version Info
%
pFooter;


