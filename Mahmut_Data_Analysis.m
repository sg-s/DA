% Mahmut_Data_Analyis.m
% 
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

% load data
load('/local-data/DA-paper/mahmut_data.mat')

redo_bootstrap = 0;

for td = 1:length(data)



 	% detrend PID with a quadratic term
	% ptrend = fit(data(td).time(:),data(td).PID(:),'Poly2'); 
	% data(td).PID = data(td).PID(:) - (ptrend(data(td).time) - mean(ptrend(data(td).time)));


	% build a simple linear model
	[K,~,filtertime] = FindBestFilter(data(td).PID(500:end),data(td).ORN(500:end),[],'filter_length=201;','min_cutoff = 0;');
	data(td).K = K;
	data(td).filtertime = filtertime*mean(diff(data(td).time));
	data(td).LinearFit = mean(data(td).ORN) + convolve(data(td).time,data(td).PID,data(td).K,data(td).filtertime);
	data(td).LinearFit(data(td).LinearFit < 0) = 0;


	figure('outerposition',[0 0 1000 400],'PaperUnits','points','PaperSize',[1000 400]); hold on
	subplot(2,1,1), hold on
	plot(data(td).time,data(td).PID,'k');
	set(gca,'XLim',[mean(data(td).time)-5 mean(data(td).time)+5])

	ylabel('PID (V)')

	subplot(2,1,2), hold on
	plot(data(td).time,data(td).ORN,'k');
	plot(data(td).time,data(td).LinearFit,'r');
	set(gca,'XLim',[mean(data(td).time)-5 mean(data(td).time)+5])
	ylabel('Firing Rate (Hz)')
	xlabel('Time (s)')
	legend ORN LinearFit
	PrettyFig;

	snapnow;
	delete(gcf);


	clear ph
	figure('outerposition',[0 0 1300 400],'PaperUnits','points','PaperSize',[1300 400]); hold on
	ph(1)=subplot(1,3,1); hold on
	ph(2)=subplot(1,3,2); hold on
	[act] = PlotDataStatistics(data,td,ph);

	subplot(1,3,3), hold on;
	plot(data(td).filtertime,data(td).K,'k','LineWidth',2)
	set(gca,'XLim',[min(data(td).filtertime) max(data(td).filtertime)],'box','on')
	xlabel('Filter Lag (s)')
	title('Filter')
	PrettyFig;

	snapnow;
	delete(gcf);


	% do gain analysis
	figure('outerposition',[0 0 1000 400],'PaperUnits','points','PaperSize',[1300 400]); hold on
	clear ph
	ph(3)=subplot(1,2,1); hold on; 	axis square
	ph(4)=subplot(1,2,2); hold on;	axis square
	s = 300; % when we start for the gain analysis
	z = length(data(td).ORN); % where we end
	step_size = 2*act/10;
	if step_size < 0.03
		step_size= 0.03;
	else
		step_size = 0.03*floor(step_size/0.03);
	end
	history_lengths = [0:step_size:2*act];
	example_history_length = history_lengths(3);

	clear x
	x.response = data(td).ORN(s:z);
	x.prediction = data(td).LinearFit(s:z);
	x.stimulus = data(td).PID(s:z);
	x.time = data(td).time(s:z);
	x.filter_length = 201;


	if redo_bootstrap
		data(td).LinearFit_p = GainAnalysis3(x,history_lengths,example_history_length,ph);
	else
		GainAnalysis3(x,history_lengths,example_history_length,ph,NaN*history_lengths);
	end
	clear x




	snapnow;
	delete(gcf);


 end
 clear td

save('/local-data/DA-paper/mahmut_data.mat','data','-append')
