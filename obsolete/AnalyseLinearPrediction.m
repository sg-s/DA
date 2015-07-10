% AnalyseLinearPrediction.m
% analyses the linear prediction on carlotta's data
function [output_data] = AnalyseLinearPrediction(filename,v,history_lengths)
if nargin < 1
	filename = '/data/random-stim/final_2012_06_05_ab3A_d2succ_randstim_100ms_60sec.mat'
end
if nargin < 2
	v = 1; % verbosity
end
if nargin < 3
	history_lengths = [30 102 150 201 300 402 501 600 801 1002];
end
% get the data
disp('Prepping data...')
if exist('mv') ~= 1
	[mv] = analysis_old_data(filename);
end
disp('DONE')

% make a figure showing the stimulus statistics
if v
	figure, subplot(1,2,1), hold on
	hist(mv(1).on_durations,15)
	xlabel('On Durations (s)','FontSize',20)

	subplot(1,2,2), hold on
	hist(mv(1).off_durations,15)
	xlabel('Off Durations (s)','FontSize',20)
	suptitle(strrep(filename,'_','-'),20)

end

output_data.ShortestPulse = min(mv(1).on_durations);
output_data.LongestPulse = max(mv(1).on_durations);

% % make a figure showing the PID
% figure, hold on
% subplot(2,1,1), hold on
% plot(t,PID,'LineWidth',1)
% plot(t,cf_PID(t),'k')
% xlabel('Time (s)')


% subplot(2,1,2), hold on
% plot(t,PID,'Color',[0.8 0.8 0.8])
% hold on
% temp = PID;
% temp(PID > 0.25 & PID < 0.75) = NaN;
% plot(t,temp,'Color',[0.4 0.4 0.4])
% temp = PID;
% temp(PID > 0.1 & PID < 0.9) = NaN;
% plot(t,temp,'r')


% translate
t=  mv(1).time;
stim = mv(1).stim;
PID = mv(1).PID;
f = mv(1).f;

% make a figure showing the data
if v
	figure, hold on
	subplot(3,1,1), hold on
	plot(t,stim,'LineWidth',2)
	set(gca,'box','on','LineWidth',2,'XLim',[15 22],'YLim',[-0.1 1.1])
	ylabel('stimulus')

	subplot(3,1,2), hold on
	plot(t,PID,'LineWidth',2)
	set(gca,'box','on','LineWidth',2,'XLim',[15 22])
	ylabel('PID (a.u.)')

	subplot(3,1,3), hold on
	plot(t,f,'LineWidth',2)
	set(gca,'box','on','LineWidth',2,'XLim',[15 22])
	xlabel('t (s)','FontSize',20)
	ylabel('f (Hz)')
	suptitle(strrep(filename,'_','-'),20)
end

% 1. now calcualte the filter and make predictions using Carloitta's code
marker_size = 10;

% translate and import variables
K1 = mv(1).K;
fitN = mv(1).NLfunction;
fp1 = mv(1).fp;
t = mv(1).time;
tfilter1 = mv(1).filtertime;
f = mv(1).f;
fs = mv(1).fs;

% NL model -- linear function
fp2 = mv(2).fp;




% plot the filter, the static nonlinearity, and the prediciton
if v
	figure, hold on,
	subplot(2,2,1), hold on
	plot(tfilter1,K1,'LineWidth',2)
	set(gca,'LineWidth',2,'FontSize',20,'box','on')
	title('Filter, PID -> f','FontSize',20)

	subplot(2,2,2), hold on
	fplot(fitN,mv(1).NLfunctionRange);
	set(gca,'LineWidth',2,'FontSize',20,'box','on')
	title('Static Nonlinearity','FontSize',20)

	subplot(2,2,3:4), hold on
	plot(t,f,'b','LineWidth',2)
	hold on
	plot(t,fp1,'r','LineWidth',2)
	plot(t,fp2,'k','LineWidth',2)
	set(gca,'LineWidth',2,'FontSize',20,'box','on')
	set(gca,'XLim',[10 15],'YLim',[0 120])
	xlabel('Time (s)')
	ylabel('Firing rate (Hz)')
	suptitle(strrep(filename,'_','-'),20)
end


% 2. now compute shat, which is the mean stimulus averaged over a window of length h at each
% time point
% we shall do this for the following history lengths:

hl = history_lengths/3; % history lengths better be divisble by 3!
shat = NaN(length(history_lengths),length(PID));
for i = 1:length(history_lengths)
	for j = (hl(i)+1):length(PID)
		shat(i,j) = mean(PID(j-hl(i):j));
	end
end


% now segment data into bottom 10% and top 10% of mean average history stimulus and plot it

% make data vectors
low_slopes = NaN*history_lengths;
low_slopes_err = NaN*history_lengths;
high_slopes = NaN*history_lengths;
high_slopes_err = NaN*history_lengths;
low_gof = NaN*history_lengths;
high_gof = NaN*history_lengths;



% prep summary figure
if v 
	summ_fig = figure; hold on;
end

for k= 1:2
	switch k
		case 1
			thisfp  = fp1; % LN model prediction
		case 2
			thisfp = fp2; % linear prediction
		
	end
	% calcualte the slopes for all points
	[fall gof] = fit(f(1:end-34),thisfp(1:end-34),'Poly1');
	all_gof = gof.rsquare;
	er = confint(fall);
	all_slopes_err=fall.p1-er(1,1);
	all_slopes = fall.p1;
	output_data.all_slopes = all_slopes;

	if v 
		figure, hold on
	end
	for i = 1:length(history_lengths)
		if v 
			subplot(2,5,i), hold on
		end
		[sorted_shat idx] = sort(shat(i,:),'ascend');
		f_low = f(idx(1:floor(length(PID)/10)));
		fp_low = thisfp(idx(1:floor(length(PID)/10)));
		[sorted_shat idx] = sort(shat(i,:),'descend');
		f_high = f(idx(1:floor(length(PID)/10)));
		fp_high = thisfp(idx(1:floor(length(PID)/10)));
		if v
			scatter(f,thisfp,marker_size,[0.7 0.7 0.7],'filled')
			scatter(f_low,fp_low,marker_size,[0.5 1 0.5],'filled')
			scatter(f_high,fp_high,marker_size,[1 0.5 0.5],'filled')
			set(gca,'LineWidth',2,'box','on','FontSize',20)
			set(gca,'XLim',[min(f)-10,max(f)+10],'YLim',[min(f)-10,max(f)+10])
			title(strcat(mat2str(history_lengths(i)),'ms'))
		end

		% remove NaN values
		f_high(isnan(fp_high)) = [];
		fp_high(isnan(fp_high)) = [];
		f_low(isnan(fp_low)) = [];
		fp_low(isnan(fp_low)) = [];

		% fit lines
		[flow gof] = fit(f_low,fp_low,'Poly1');
		low_gof(i) = gof.rsquare;
		er = confint(flow);
		low_slopes_err(i)=flow.p1-er(1,1);
		low_slopes(i) = flow.p1;

		[fhigh gof] = fit(f_high,fp_high,'Poly1');
		high_gof(i) =  gof.rsquare;
		er = confint(fhigh);
		high_slopes_err(i)=fhigh.p1-er(1,1);
		high_slopes(i) = fhigh.p1;



		% plot the lines
		if v
			plot(min(f):max(f),fall(min(f):max(f)),'Color',[0.5 0.5 0.5],'LineWidth',2)
			plot(min(f_low):max(f_low),flow(min(f_low):max(f_low)),'g','LineWidth',2)
			plot(min(f_high):max(f_high),fhigh(min(f_high):max(f_high)),'r','LineWidth',2)
			suptitle(mat2str(k))
		end
		
		

	end
	

	% plot to summary figure
	if v
		figure(summ_fig), hold on
		subplot(2,2,k*2-1), hold on
		plot(history_lengths,all_slopes*ones(1,length(history_lengths)),'k','LineWidth',2), hold on
		errorbar(history_lengths,low_slopes,low_slopes_err,'b','LineWidth',2), hold on
		errorbar(history_lengths,high_slopes,high_slopes_err,'r','LineWidth',2)
		set(gca,'LineWidth',2,'FontSize',20,'box','on')
		xlabel('History Length (ms)','FontSize',20)
		ylabel('Slope','FontSize',20)

		subplot(2,2,k*2), hold on
		plot(history_lengths,low_gof,'b','LineWidth',2),hold on
		plot(history_lengths,high_gof,'r','LineWidth',2)
		xlabel('History Length (ms)','FontSize',20)
		ylabel('Goodness of Fit, rsquare','FontSize',20)
		set(gca,'LineWidth',2,'FontSize',20,'box','on')


	end

	output_data(k).high_slopes = high_slopes;
	output_data(k).high_slopes_err = high_slopes_err;
	output_data(k).low_slopes = low_slopes;
	output_data(k).low_slopes_err = low_slopes_err;
	output_data(k).low_gof = low_gof;
	output_data(k).high_gof = high_gof;


end
if v 
	suptitle(strrep(filename,'_','-'),20)
end


