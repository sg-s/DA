
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


%% Version Info
%
pFooter;


