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

function [] = GainAnalysis3(x,history_lengths,example_history_length,plothere)

% set defaults
marker_size=10;
marker_size2=24;
font_size=20;
plotid=[1 2];

if nargin < 4
	figure, hold on; 
	plothere(1) = subplot(2,1,1); hold on;
	plothere(2) = subplot(2,1,2); hold on;
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
shat = NaN(length(hl),length(stimulus));
for i = 1:length(hl)
	shat(i,:) = filter(ones(1,hl(i))/hl(i),1,stimulus);
	shat(i,1:hl(i)) = NaN;
end



% make data vectors
low_slopes = NaN*history_lengths;
low_slopes_err = NaN*history_lengths;
high_slopes = NaN*history_lengths;
high_slopes_err = NaN*history_lengths;
low_gof = NaN*history_lengths;
high_gof = NaN*history_lengths;



% calculate the slopes for all points
[fall gof] = fit(fp(filter_length+2:end-filter_length),f(filter_length+2:end-filter_length),'Poly1');
all_gof = gof.rsquare;
er = confint(fall);
all_slopes_err=fall.p1-er(1,1);
all_slopes = fall.p1;
output_data.all_slopes = all_slopes;



for i = 1:length(history_lengths)

	% find lowest and highest 10%
	this_shat = shat(i,:);
	this_shat(1:hl(i)) = Inf; % the initial segment where we can't estimate shat is excluded
	[sorted_shat idx] = sort(this_shat,'ascend');
	f_low = f(idx(1:floor(length(stimulus)/10)));
	fp_low = fp(idx(1:floor(length(stimulus)/10)));


	this_shat(1:hl(i)) = -Inf;
	[sorted_shat idx] = sort(this_shat,'descend');
	f_high = f(idx(1:floor(length(stimulus)/10)));
	fp_high = fp(idx(1:floor(length(stimulus)/10)));

	% remove NaN values
	f_high(isnan(fp_high)) = [];
	fp_high(isnan(fp_high)) = [];
	f_low(isnan(fp_low)) = [];
	fp_low(isnan(fp_low)) = [];

	% fit lines
	[flow gof] = fit(fp_low,f_low,'Poly1');
	low_gof(i) = gof.rsquare;
	er = confint(flow);
	low_slopes_err(i)=flow.p1-er(1,1);
	low_slopes(i) = flow.p1;

	[fhigh gof] = fit(fp_high,f_high,'Poly1');
	high_gof(i) =  gof.rsquare;
	er = confint(fhigh);
	high_slopes_err(i)=fhigh.p1-er(1,1);
	high_slopes(i) = fhigh.p1;

	if history_lengths(i) == example_history_length
		disp(example_history_length)

		% plot the stimulus and the smoothed stimulus
		plot(plothere(1),t,stimulus,'k','LineWidth',2), hold on
		plot(plothere(1),t,shat(i,:),'Color',[0.9 0.9 0.9],'LineWidth',4)

		% plot the response and the prediction 
		plot(plothere(2),t,f,'k','LineWidth',2), hold on
		plot(plothere(2),t,fp,'r','LineWidth',2), hold on

		
		% indicate regions of lowest and highest 10%
		tp = floor(length(stimulus)/10);
		scatter(plothere(1),t(idx(1:tp)),shat(i,idx(1:tp)),'r','fill')
		scatter(plothere(1),t(idx(end-tp:end)),shat(i,idx(end-tp:end)),'g','fill')
		% scatter(plothere(2),t(idx(1:tp)),f(idx(1:tp)),'r','fill')
		scatter(plothere(2),t(idx(1:tp)),fp(idx(1:tp)),'r','fill')
		scatter(plothere(2),t(idx(end-tp:end)),fp(idx(end-tp:end)),'g','fill')
		% scatter(plothere(2),t(idx(end-tp:end)),f(idx(end-tp:end)),'g','fill')

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

	end	

end


% plot to summary figure
plot(plothere(4),history_lengths,all_slopes*ones(1,length(history_lengths)),'k','LineWidth',2), hold on
errorbar(plothere(4),history_lengths,low_slopes,low_slopes_err,'g','LineWidth',2), hold on
errorbar(plothere(4),history_lengths,high_slopes,high_slopes_err,'r','LineWidth',2)
	
% bootstrap slopes
% [low_slopes, high_slopes, p] = BootStrapErrorBars(x,history_lengths,0.1);
% plot(history_lengths,low_slopes.data,'g','LineWidth',2), hold on
% plot(history_lengths,high_slopes.data,'r','LineWidth',2)
% sig = p<(0.05/length(p)); % these points are significant, Bonferroni corrected
% scatter(history_lengths(sig),low_slopes.data(sig),1256,'g.')
% scatter(history_lengths(sig),high_slopes.data(sig),1256,'r.')

return

	set(gca,'LineWidth',2,'FontSize',20,'box','on','XLim',[0 max(history_lengths)])
	xlabel('History Length (ms)','FontSize',20)
	ylabel('Slope data/prediction (gain)','FontSize',20)

	if ismember(3,plotid)
		subplot(1,2,2), hold on
		plot(history_lengths,low_gof,'g','LineWidth',2),hold on
		plot(history_lengths,high_gof,'r','LineWidth',2)
		xlabel('History Length (ms)','FontSize',20)
		ylabel('Goodness of Fit, rsquare','FontSize',20)
		set(gca,'LineWidth',2,'FontSize',20,'box','on','YLim',[0 1.1],'XLim',[0 max(history_lengths)])
	end


output_data.high_slopes = high_slopes;
output_data.high_slopes_err = high_slopes_err;
output_data.low_slopes = low_slopes;
output_data.low_slopes_err = low_slopes_err;
output_data.low_gof = low_gof;
output_data.high_gof = high_gof;