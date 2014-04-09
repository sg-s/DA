function [output_data] = GainAnalysis(f,fp,PID,shat,history_lengths,hl,filter_length,marker_size,marker_size2,font_size,plotid,plothere)

if nargin < 12
	plothere = [];
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
	f_low = f(idx(1:floor(length(PID)/10)));
	fp_low = fp(idx(1:floor(length(PID)/10)));

	this_shat(1:hl(i)) = -Inf;
	[sorted_shat idx] = sort(this_shat,'descend');
	f_high = f(idx(1:floor(length(PID)/10)));
	fp_high = fp(idx(1:floor(length(PID)/10)));

	if ismember(1,plotid)
		plot(fp,f,'.','MarkerSize',marker_size,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
		plot(fp_low,f_low,'.','MarkerSize',marker_size,'MarkerFaceColor',[0.5 1 0.5],'MarkerEdgeColor',[0.5 1 0.5])
		plot(fp_high,f_high,'.','MarkerSize',marker_size,'MarkerFaceColor',[1 0.5 0.5],'MarkerEdgeColor',[1 0.5 0.5])
		set(gca,'LineWidth',2,'box','on','FontSize',font_size)
		set(gca,'XLim',[min(f)-10,max(f)+10],'YLim',[min(f)-10,max(f)+10])
		axis square
		title(strcat(mat2str(history_lengths(i)),'ms'))

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
		plot(min(f):max(f),fall(min(f):max(f)),'Color',[0.5 0.5 0.5],'LineWidth',2)
		plot(min(f_low):max(f_low),flow(min(f_low):max(f_low)),'g','LineWidth',2)
		plot(min(f_high):max(f_high),fhigh(min(f_high):max(f_high)),'r','LineWidth',2)
	end

end



if ismember(2,plotid)
	% plot to summary figure
	figure('outerposition',[10 10 1000 500],'PaperUnits','points','PaperSize',[1000 700]); hold on
	subplot(1,2,1), hold on

	
	plot(history_lengths,all_slopes*ones(1,length(history_lengths)),'k','LineWidth',2), hold on
	errorbar(history_lengths,low_slopes,low_slopes_err,'g','LineWidth',2), hold on
	errorbar(history_lengths,high_slopes,high_slopes_err,'r','LineWidth',2)
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