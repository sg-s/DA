% GainAnalysis3.m
% a better, improved version of GainAnalysis.m
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
% performs a gain analysis and makes some plots

function [p,low_slopes,high_slopes,low_gof,high_gof] = GainAnalysis3(x,history_lengths,example_history_length,plothere,p)


% set defaults
marker_size=10;
marker_size2=24;
font_size=20;
plotid=[1 2 3 4];
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

	% find highest 10%
	this_shat = shat(i,:);

	% assign undesirables to -Inf
	this_shat(1:hl(i)) = -Inf;
	this_shat(isnan(this_shat)) = -Inf;
	[~, idx] = sort(this_shat,'descend');
	f_high = f(idx(1:n));
	fp_high = fp(idx(1:n));

	% remove NaN values
	censor_these = find(isnan(fp_high) + isnan(fp_low) + isnan(f_low) + isnan(f_high));
	f_high(censor_these) = [];
	fp_high(censor_these) = [];
	f_low(censor_these) = [];
	fp_low(censor_these) = [];


	% censor times when f is 0?
	f_low(f==0) = [];
	f_high(f==0) = [];
	f(f==0) = [];

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
			plot(plothere(2),t,f,'k','LineWidth',2), hold on
			plot(plothere(2),t,fp,'r','LineWidth',2), hold on

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
			title(plothere(3),strcat(mat2str(history_lengths(i)),'s'))

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
	[low_slopes2, high_slopes2, p] = BootStrapErrorBars(x,history_lengths,0.1);
end

if length(plothere) == 4
	plot(plothere(4),history_lengths,low_slopes2.data,'g','LineWidth',2), hold on
	plot(plothere(4),history_lengths,high_slopes2.data,'r','LineWidth',2)
	p = p*length(p); % Bonferroni correction
	sig = p<0.05; % these points are significant,
	scatter(plothere(4),history_lengths(sig),low_slopes2.data(sig),1256,'g.')
	scatter(plothere(4),history_lengths(sig),high_slopes2.data(sig),1256,'r.')

	set(plothere(4),'LineWidth',2,'FontSize',20,'box','on','XLim',[0 max(history_lengths)])
	xlabel(plothere(4),'History Length (s)','FontSize',20)
	ylabel(plothere(4),'Slope data/prediction (gain)','FontSize',20)
end


% clean up slopes to send out
low_slopes = low_slopes2.data;
high_slopes = high_slopes2.data;