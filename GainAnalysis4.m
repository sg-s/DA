% GainAnalysis4.m
% a better, improved version of GainAnalysis.m
% GainAnalysis4 differs from GainAnalysis3 in that p_values are calculated for each cloud, instead of pair-wise
%
% x is a structure with the following fields:
% x.response -- the actual data
% x.prediction -- the prediction you want to compare to
% x.time -- time vector
% x.stimulus -- the stimulus vector 
%
% and you also need to specify
%
% history_lengths (a vector of history lengths)
% example_history_length (one number)
% plothere is a 4x1 vector of axis handles if you want to plots to appear somewhere special. the first two are used for fig 1 
%
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [p,low_slopes,high_slopes,low_gof,high_gof,example_plot,extra_variables] = GainAnalysis4(x,history_lengths,example_history_length,plothere,p)


% set defaults
example_plot = [];
marker_size=10;
marker_size2=24;
font_size=20;
plotid=[1 2 3 4]; %specify which plot you want
debug = 0;

switch nargin 
case 1
	error('Need to specify which history lengths to analyse')
case 2
	% don't plot anything
	plotid = [];
	plothere = [];
	example_history_length = [];
case 3
case 4
case 5
	% split p into p_low and p_high
	p_low = p(:,1);
	p_high = p(:,2);

end


if isempty(plothere) && ~isempty(plotid)
	if debug
		figure, hold on; 
		plothere(1) = subplot(2,1,1); hold on;
		plothere(2) = subplot(2,1,2); hold on;
	end
	figure, hold on;
	plothere(3) = gca;
	figure; hold on;
	plothere(4) = gca;
end


% unpack data
f = x.response(:);
fp = x.prediction(:);
stimulus = x.stimulus(:);
t = x.time(:);
filter_length = x.filter_length;


% figure out sampling rate
dt = mean(diff(t));
hl = round(history_lengths/dt);



% compute shat(the smoothed stimulus)
shat = ComputeSmoothedStimulus(stimulus,hl);

% make data vectors
low_slopes = NaN*history_lengths;
low_slopes_err = NaN*history_lengths;
high_slopes = NaN*history_lengths;
high_slopes_err = NaN*history_lengths;
low_gof = NaN*history_lengths;
high_gof = NaN*history_lengths;

% some extra data vectors
low_min = NaN*history_lengths;
low_max = NaN*history_lengths;
high_min = NaN*history_lengths;
high_max = NaN*history_lengths;


% calculate the slopes for all points
[fall, gof] = fit(fp(filter_length+2:end-filter_length),f(filter_length+2:end-filter_length),'Poly1');
all_gof = gof.rsquare;
er = confint(fall);
all_slopes_err=fall.p1-er(1,1);
all_slopes = fall.p1;
output_data.all_slopes = all_slopes;

n = floor(sum(~isnan(stimulus))/10);


for i = 1:length(history_lengths)

	f = x.response(:);

	% find lowest 10%
	this_shat = shat(i,:);

	% assign undesirables to Inf
	this_shat(1:hl(i)) = Inf; % the initial segment where we can't estimate shat is excluded
	this_shat(isnan(this_shat)) = Inf;
	[~, idx] = sort(this_shat,'ascend');
	f_low = f(idx(1:n));
	fp_low = fp(idx(1:n));
	t_low = t(idx(1:n));
	
		

	% find highest 10%
	this_shat = shat(i,:);

	% assign undesirables to -Inf
	this_shat(1:hl(i)) = -Inf;
	this_shat(isnan(this_shat)) = -Inf;
	[~, idx] = sort(this_shat,'descend');
	f_high = f(idx(1:n));
	fp_high = fp(idx(1:n));
	t_high = t(idx(1:n));

	% remove NaN values
	censor_these = find(isnan(fp_high) + isnan(fp_low) + isnan(f_low) + isnan(f_high));
	f_high(censor_these) = [];
	fp_high(censor_these) = [];
	f_low(censor_these) = [];
	fp_low(censor_these) = [];
	t_low(censor_these) = [];
	t_high(censor_these) = [];

	% determine range of this subset of data
	low_min(i) = min(f_low);
	low_max(i) = max(f_low);
	high_min(i) = min(f_high);
	high_max(i) = max(f_high);

	% fit lines
	[flow, gof] = fit(fp_low,f_low,'Poly1');
	low_gof(i) = gof.rsquare;
	er = confint(flow);
	low_slopes_err(i)=flow.p1-er(1,1);
	low_slopes(i) = flow.p1;

	[fhigh, gof] = fit(fp_high,f_high,'Poly1');
	high_gof(i) =  gof.rsquare;
	er = confint(fhigh);
	high_slopes_err(i)=fhigh.p1-er(1,1);
	high_slopes(i) = fhigh.p1;

	

	if history_lengths(i) == example_history_length


		if debug || plothere(1)
			% % plot the stimulus and the smoothed stimulus
			plot(plothere(1),t,stimulus,'k','LineWidth',2), hold on
			plot(plothere(1),t,shat(i,:),'Color',[0.9 0.9 0.9],'LineWidth',4)

			% plot the response and the prediction 
			plot(plothere(2),t,f,'k','LineWidth',1), hold on
			plot(plothere(2),t,fp,'r','LineWidth',1), hold on

			% highlight the sections we use for the analysis
			plot(plothere(2),t_high,f_high,'k.','LineWidth',2), hold on
			plot(plothere(2),t_high,fp_high,'r.','LineWidth',2), hold on
			plot(plothere(2),t_low,f_low,'k.','LineWidth',2), hold on
			plot(plothere(2),t_low,fp_low,'r.','LineWidth',2), hold on

			% save these for later
			example_plot.t_low = t_low;
			example_plot.f_low = f_low;
			example_plot.fp_low = fp_low;
			example_plot.t_high= t_high;
			example_plot.f_high= f_high;
			example_plot.fp_high = fp_high;

			xlabel(plothere(2),'Time (s)','FontSize',20)
			xlabel(plothere(1),'Time (s)','FontSize',20)

			ylabel(plothere(2),'Firing rate (Hz)','FontSize',20)
			ylabel(plothere(1),'Stimulus (a.u.)','FontSize',20)
			
			% indicate regions of lowest and highest 10%
			tp = floor(length(stimulus)/10);
			plot(plothere(1),t(idx(1:tp)),shat(i,idx(1:tp)),'r.')
			plot(plothere(1),t(idx(end-tp:end)),shat(i,idx(end-tp:end)),'g.')


			set(plothere(1),'LineWidth',2,'FontSize',20)
			set(plothere(2),'LineWidth',2,'FontSize',20)

		end

		if plothere(3)
			% plot these on the scatter plot
			plot(plothere(3),fp,f,'.','MarkerSize',marker_size,'MarkerFaceColor',[0.9 0.9 0.9],'MarkerEdgeColor',[0.9 0.9 0.9])
			plot(plothere(3),fp_low,f_low,'.','MarkerSize',marker_size,'MarkerFaceColor',[0.5 1 0.5],'MarkerEdgeColor',[0.5 1 0.5])
			plot(plothere(3),fp_high,f_high,'.','MarkerSize',marker_size,'MarkerFaceColor',[1 0.5 0.5],'MarkerEdgeColor',[1 0.5 0.5])
			set(plothere(3),'LineWidth',2,'box','on','FontSize',font_size)
			set(plothere(3),'XLim',[min(f)-1,max(f)+1],'YLim',[min(f)-1,max(f)+1])

			titlestr = strcat('\tau_h=',mat2str(history_lengths(i)),'s  gof:');
			titlestr = strcat(titlestr, '{\color{green}',oval(low_gof(i),2),'},');
			titlestr = strcat(titlestr, '{\color{red}',oval(high_gof(i),2),'}');

			title(plothere(3),titlestr);

			% plot the best fit lines
			plot(plothere(3),[min(fp) max(fp)],fall([min(fp) max(fp)]),'Color',[0.5 0.5 0.5],'LineWidth',3)
			plot(plothere(3),[min(fp_low) max(fp_low)],flow([min(fp_low) max(fp_low)]),'g','LineWidth',3)
			plot(plothere(3),[min(fp_high) max(fp_high)],fhigh([min(fp_high) max(fp_high)]),'r','LineWidth',3)

			xlabel(plothere(3),'Prediction')
			ylabel(plothere(3),'Actual neuron response')
		end

	end	

end

if length(plothere) == 4
	% plot to summary figure
	plot(plothere(4),history_lengths,all_slopes*ones(1,length(history_lengths)),'k','LineWidth',2), hold on
end

	
% bootstrap slopes
if nargin == 5
	low_slopes2.data = low_slopes;
	high_slopes2.data = high_slopes;

else
	[low_slopes2, high_slopes2] = BootStrapErrorBars(x,history_lengths,0.1);

	% calculate p-values for each cloud 
	p_low  = ones(length(history_lengths),1);
	p_high = ones(length(history_lengths),1);
	for i = 1:length(history_lengths)
		a = abs(low_slopes2.bootstrap(:,i) - all_slopes);
		a0 = abs(low_slopes2.data(i) - all_slopes);
		p_low(i) = sum(a>a0)/length(low_slopes2.bootstrap);

		a = abs(high_slopes2.bootstrap(:,i) - all_slopes);
		a0 = abs(high_slopes2.data(i) - all_slopes);
		p_high(i) = sum(a>a0)/length(high_slopes2.bootstrap);
	end

end

if length(plothere) == 4
	plot(plothere(4),history_lengths,low_slopes2.data,'g','LineWidth',2), hold on
	plot(plothere(4),history_lengths,high_slopes2.data,'r','LineWidth',2)
	%p_low = p_low*length(p_low); % Bonferroni correction
	%p_high = p_high*length(p_high); % Bonferroni correction
	sig_low = p_low<0.05; % these points are significant,
	sig_high = p_high<0.05; % these points are significant,
	scatter(plothere(4),history_lengths(sig_low),low_slopes2.data(sig_low),1256,'g.')
	scatter(plothere(4),history_lengths(sig_high),high_slopes2.data(sig_high),1256,'r.')

	set(plothere(4),'LineWidth',2,'FontSize',20,'box','on','XLim',[0 max(history_lengths)])
	xlabel(plothere(4),'History Length (s)','FontSize',20)
	ylabel(plothere(4),'Slope data/prediction (gain)','FontSize',20)
end


% clean up slopes to send out
low_slopes = low_slopes2.data;
high_slopes = high_slopes2.data;

% package p for backwards compatibility
p = [p_low p_high];

% package some extra variables that may be requested
extra_variables.data_min = min(x.response)*ones(length(history_lengths),1);
extra_variables.data_max = max(x.response)*ones(length(history_lengths),1);
extra_variables.low_min = low_min;
extra_variables.low_max = low_max;
extra_variables.high_min = high_min;
extra_variables.high_max = high_max;

extra_variables.all_slopes = fall.p1*ones(length(history_lengths),1);