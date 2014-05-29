% GainAnalysis2.m
% a better, improved version of GainAnalysis.m
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
% performs a gain analysis and makes some plots
% x is a structure with the following fields:
% x.data -- the actual data
% x.prediction -- the prediction you want to compare to
% x.time -- time vector
% x.stimulus -- the stimulus vector 
%

function [output_data] = GainAnalysis2(x,history_lengths,filter_length,varargin)

% set defaults
marker_size=10;
marker_size2=24;
font_size=20;
plotid=[1 2];
plothere=[];

if iseven(nargin-3)
    for i = 1:length(varargin)
        temp = varargin{i};
        if ischar(temp)
            eval(strcat(temp,'=varargin{i+1};'));
        end
    end
    
else
    error('Inputs need to be name value pairs')
end
	

f = x.data;
fp = x.prediction;
stimulus = x.stimulus;
t = x.time;


% figure out sampling rate
dt = mean(diff(t));
hl = round(history_lengths/dt);


% compute shat
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

f = f(:);
fp = fp(:);
[fall gof] = fit(fp(filter_length+2:end),f(filter_length+2:end),'Poly1');
all_gof = gof.rsquare;
er = confint(fall);
all_slopes_err=fall.p1-er(1,1);
all_slopes = fall.p1;
output_data.all_slopes = all_slopes;
shuffled_low = zeros(1,length(history_lengths));
shuffled_high = zeros(1,length(history_lengths));


if ismember(1,plotid) && isempty(plothere)
	f3= figure('outerposition',[0 0 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
end
for i = 1:length(history_lengths)

	if ismember(1,plotid) && isempty(plothere)
		figure(f3), hold on
		autoplot(length(hl),i); hold on
	elseif ismember(1,plotid) 
		axis(plothere); hold on
	end

	this_shat = shat(i,:);
	this_shat(1:hl(i)) = Inf; % the initial segment where we can't estimate shat is excluded
	[sorted_shat idx] = sort(this_shat,'ascend');
	f_low = f(idx(1:floor(length(stimulus)/10)));
	fp_low = fp(idx(1:floor(length(stimulus)/10)));


	this_shat(1:hl(i)) = -Inf;
	[sorted_shat idx] = sort(this_shat,'descend');
	f_high = f(idx(1:floor(length(stimulus)/10)));
	fp_high = fp(idx(1:floor(length(stimulus)/10)));

	if ismember(1,plotid)
		plot(fp,f,'.','MarkerSize',marker_size,'MarkerFaceColor',[0.9 0.9 0.9],'MarkerEdgeColor',[0.9 0.9 0.9])
		plot(fp_low,f_low,'.','MarkerSize',marker_size,'MarkerFaceColor',[0.5 1 0.5],'MarkerEdgeColor',[0.5 1 0.5])
		plot(fp_high,f_high,'.','MarkerSize',marker_size,'MarkerFaceColor',[1 0.5 0.5],'MarkerEdgeColor',[1 0.5 0.5])
		set(gca,'LineWidth',2,'box','on','FontSize',font_size)
		set(gca,'XLim',[min(f)-1,max(f)+1],'YLim',[min(f)-1,max(f)+1])
		axis square
		title(strcat(mat2str(history_lengths(i)),'s'))

		if i == 6
			xlabel('Linear Prediction (Hz)','FontSize',16)
			ylabel('Neuron Response (Hz)','FontSize',16)
		end
	end

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

	if ismember(1,plotid)

		% plot the best fit lines
		plot([min(fp) max(fp)],fall([min(fp) max(fp)]),'Color',[0.5 0.5 0.5],'LineWidth',3)
		plot([min(fp_low) max(fp_low)],flow([min(fp_low) max(fp_low)]),'g','LineWidth',3)
		plot([min(fp_high) max(fp_high)],fhigh([min(fp_high) max(fp_high)]),'r','LineWidth',3)
	end

	
end



if ismember(2,plotid)
	% plot to summary figure
	subplot(plothere)
	plot(history_lengths,all_slopes*ones(1,length(history_lengths)),'k','LineWidth',2), hold on
	errorbar(history_lengths,low_slopes,low_slopes_err,'g','LineWidth',2), hold on
	errorbar(history_lengths,high_slopes,high_slopes_err,'r','LineWidth',2)

	% bootstrap slopes
	%[low_slopes_min, low_slopes_max, high_slopes_min, high_slopes_max] = BootStrapErrorBars(shat,f,fp,hl,fraction);

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
end

output_data.high_slopes = high_slopes;
output_data.high_slopes_err = high_slopes_err;
output_data.low_slopes = low_slopes;
output_data.low_slopes_err = low_slopes_err;
output_data.low_gof = low_gof;
output_data.high_gof = high_gof;