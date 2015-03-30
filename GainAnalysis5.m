% GainAnalysis5.m
% a better, improved version of GainAnalysis4.m
% see: https://github.com/sg-s/DA/issues/182
%
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [p,low_slopes,high_slopes,low_gof,high_gof,example_plot,extra_variables] = GainAnalysis5(x,history_lengths,example_history_length,plothere,p)


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
	p_low = p(1,:);
	p_high = p(2,:);
case 6
	% split p into p_low and p_high
	p_low = p(1,:);
	p_high = p(2,:);

end


if isempty(plothere) 
	if debug
		figure, hold on; 
		plothere(1) = subplot(2,1,1); hold on;
		plothere(2) = subplot(2,1,2); hold on;
		figure, hold on;
		plothere(3) = gca;
		figure; hold on;
		plothere(4) = gca;
	end
	
end


% unpack data
f = x.response(:);
fp = x.prediction(:);
stimulus = x.stimulus(:);
t = x.time(:);
frac = x.frac;


% figure out sampling rate
dt = mean(diff(t));
hl = round(history_lengths/dt);

%         ######  ##     ##  #######   #######  ######## ##     ## 
%        ##    ## ###   ### ##     ## ##     ##    ##    ##     ## 
%        ##       #### #### ##     ## ##     ##    ##    ##     ## 
%         ######  ## ### ## ##     ## ##     ##    ##    ######### 
%              ## ##     ## ##     ## ##     ##    ##    ##     ## 
%        ##    ## ##     ## ##     ## ##     ##    ##    ##     ## 
%         ######  ##     ##  #######   #######     ##    ##     ## 


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
[coeff,score,latent] = pca([fp(~(isnan(fp) | isnan(f))) f(~(isnan(fp) | isnan(f)))]);
all_slopes = coeff(2,1)/coeff(1,1);





n = floor(sum(~isnan(f))*frac);

for i = 1:length(history_lengths)

	this_shat = shat(i,:);


	% find times when smoothed stim (this_shat) is in lowest 10%
	this_shat(1:hl(i)) = Inf; % the initial segment where we can't estimate shat is excluded
	this_shat(isnan(this_shat)) = Inf;
	this_shat(isnan(f)) = Inf;
	[~, t_low] = sort(this_shat,'ascend');
	t_low = t_low(1:n); % this is an index
	f_low = f(t_low);
	fp_low = fp(t_low);
	s_low = this_shat(t_low);
	t_low = t(t_low); % t_low is now a time. 
	 
	this_shat = shat(i,:);
	% find times when smoothed stim is highest 10%
	this_shat(1:hl(i)) = -Inf;
	this_shat(isinf(this_shat)) = -Inf;
	this_shat(isnan(f)) = -Inf;
	[~, t_high] = sort(this_shat,'descend');
	t_high = t_high(1:n);
	f_high = f(t_high);
	fp_high = fp(t_high);
	s_high = this_shat(t_high);
	t_high  = t(t_high);


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

	% use PCA to get slopes of clouds of points
	[coeff_low,score_low,latent] = pca([fp_low f_low]);
	low_slopes(i) = coeff_low(2,1)/coeff_low(1,1);
	low_gof(i) = latent(1)/sum(latent);



	[coeff_high,score_high,latent] = pca([fp_high f_high]);
	high_slopes(i) = coeff_high(2,1)/coeff_high(1,1);
	high_gof(i) = latent(1)/sum(latent);



	if history_lengths(i) == example_history_length
		dothis = 0;
		if ~isempty(plothere)
			temp = plothere(1);
			temp = whos('temp');
			temp = temp.class;
			dothis=0;
			if strcmp(temp,'matlab.graphics.axis.Axes')
				dothis=1;
			elseif strcmp(temp,'double')
				dothis = plothere(1);
			end
		
			if debug || dothis
				% % plot the stimulus and the smoothed stimulus
				plot(plothere(1),t,stimulus,'k','LineWidth',2), hold on
				plot(plothere(1),t,shat(i,:),'Color',[0.9 0.9 0.9],'LineWidth',4)

				% plot the response and the prediction 
				plot(plothere(2),t,f,'k','LineWidth',1), hold on
				plot(plothere(2),t,fp,'r','LineWidth',1), hold on

				plot(plothere(1),t_high,s_high,'r.')
				plot(plothere(1),t_low,s_low,'g.')

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

				ylabel(plothere(2),'Firing rate (Hz)')
				ylabel(plothere(1),'Stimulus (a.u.)')
				
				% indicate regions of lowest and highest 10%
				tp = floor(frac*length(stimulus));


			end

			temp = plothere(3);
			temp = whos('temp');
			temp = temp.class;

			dothis=0;
			if strcmp(temp,'matlab.graphics.axis.Axes')
				dothis=1;
			elseif strcmp(temp,'double')
				dothis = plothere(3);
			end
		end
	
	
		if dothis


			hold(plothere(3),'on')

			% plot these on the scatter plot
			ss = floor(length(fp_low)/1000); % plot only a 1000 points


			plot(plothere(3),fp(1:ss:end),f(1:ss:end),'.','MarkerSize',marker_size,'MarkerFaceColor',[0.9 0.9 0.9],'MarkerEdgeColor',[0.9 0.9 0.9])
			plot(plothere(3),fp_low(1:ss:end),f_low(1:ss:end),'.','MarkerSize',marker_size,'MarkerFaceColor',[0.5 1 0.5],'MarkerEdgeColor',[0.5 1 0.5])
			plot(plothere(3),fp_high(1:ss:end),f_high(1:ss:end),'.','MarkerSize',marker_size,'MarkerFaceColor',[1 0.5 0.5],'MarkerEdgeColor',[1 0.5 0.5])
			set(plothere(3),'box','on')
			set(plothere(3),'XLim',[min(f)-1,max(f)+1],'YLim',[min(f)-1,max(f)+1])

			titlestr = strcat('\tau_h=',mat2str(history_lengths(i)),'s  gof:');
			titlestr = strcat(titlestr, '{\color{green}',oval(low_gof(i),2),'},');
			titlestr = strcat(titlestr, '{\color{red}',oval(high_gof(i),2),'}');

			title(plothere(3),titlestr);

			% plot the best fit lines
			xx = coeff_low(1,1).*score_low(:,1); xx = xx+mean(fp_low); 
			[xx,idx] = sort(xx);
			y = coeff_low(2,1).*score_low(:,1); y = y+mean(f_low);
			y = y(idx);
			plot(plothere(3),xx,y,'g')

			xx = coeff_high(1,1).*score_high(:,1); xx = xx+mean(fp_high); 
			[xx,idx] = sort(xx);
			y = coeff_high(2,1).*score_high(:,1); y = y+mean(f_high);
			y = y(idx);
			plot(plothere(3),xx,y,'r')

			xlabel(plothere(3),'Prediction')
			ylabel(plothere(3),'Actual neuron response')
		end
	end	
end

if length(plothere) == 4
	% plot to summary figure
	plot(plothere(4),history_lengths,all_slopes*ones(1,length(history_lengths)),'k'), hold on
end


%    ########   #######   #######  ########  ######  ######## ########     ###    ########  
%    ##     ## ##     ## ##     ##    ##    ##    ##    ##    ##     ##   ## ##   ##     ## 
%    ##     ## ##     ## ##     ##    ##    ##          ##    ##     ##  ##   ##  ##     ## 
%    ########  ##     ## ##     ##    ##     ######     ##    ########  ##     ## ########  
%    ##     ## ##     ## ##     ##    ##          ##    ##    ##   ##   ######### ##        
%    ##     ## ##     ## ##     ##    ##    ##    ##    ##    ##    ##  ##     ## ##        
%    ########   #######   #######     ##     ######     ##    ##     ## ##     ## ##        


% bootstrap slopes
if nargin > 4 % we specify the p-values
	low_slopes2.data = low_slopes;
	high_slopes2.data = high_slopes;
	
else
	[low_slopes2, high_slopes2,p] = BootStrapErrorBars5(x,history_lengths,frac,low_slopes,high_slopes);

	p_low = p; p_high = p;

end

% throw out all invalid data
% v = (validity.H < 0.1 & validity.R_low > .5 & validity.R_high > .5 & validity.low_gof > .8 & validity.high_gof > .8);
% low_slopes2.data(~v) = NaN;
% high_slopes2.data(~v) = NaN;

	


% this is where we plot the history length plot
if length(plothere) == 4
	hold(plothere(4),'on')

	% plot line to indicate the location of the example history plot
	yy=([.8*min(high_slopes2.data) 1.2*max(low_slopes2.data)]);
	if ~isempty(example_history_length)
		plot(plothere(4),[example_history_length example_history_length],yy,'k-.')
	end

	plot(plothere(4),history_lengths,low_slopes2.data,'g'), hold on
	plot(plothere(4),history_lengths,high_slopes2.data,'r')

	% we now perform a Holm-Bonferroni correction 
	% (see: https://en.wikipedia.org/wiki/Holm%E2%80%93Bonferroni_method)

	if any(isnan(p_low))
		% this means we should skip this...
		sig_low = 0*history_lengths;
		sig_low = logical(sig_low);
		sig_high = sig_low;
	else

		% order p values from low to high
		[p_low,idx] = sort(p_low);
		history_lengths = history_lengths(idx);
		A = .05; % alpha, significance value

		% find the minimum k such that p(k) > A/(m+1-k)
		m = length(history_lengths);
		k_cut = [];
		for k = 1:m
			if p_low(k) > A/(m + 1 - k)
				k_cut = k;
				break
			end
		end
		
		% accept k-1 hypotheses
		if isempty(k_cut)
			% everything is signficant 
			sig_low = 1+0*history_lengths;
			sig_low = logical(sig_low);
		else
			sig_low = 0*history_lengths;
			sig_low(1:k_cut) = 1; 
			sig_low = logical(sig_low);
		end
		sig_high = sig_low;
		p_high = p_low;
		
	end

	% resort the history lengths into the original order
	[history_lengths,idx] = sort(history_lengths);
	sig_low = sig_low(idx);
	sig_high = sig_high(idx);
	p_low = p_low(idx);
	p_high = p_high(idx);

	scatter(plothere(4),history_lengths(sig_low),low_slopes2.data(sig_low),1256,'g.')
	scatter(plothere(4),history_lengths(sig_high),high_slopes2.data(sig_high),1256,'r.')

	

	set(plothere(4),'box','on','XLim',[0 max(history_lengths)])
	xlabel(plothere(4),'History Length (s)','FontSize',20)
	ylabel(plothere(4),'Slope data/prediction (gain)','FontSize',20)
end



% package p for backwards compatibility
p = [p_low; p_high];

% package some extra variables that may be requested
extra_variables.data_min = min(x.response)*ones(length(history_lengths),1);
extra_variables.data_max = max(x.response)*ones(length(history_lengths),1);
extra_variables.low_min = low_min;
extra_variables.low_max = low_max;
extra_variables.high_min = high_min;
extra_variables.high_max = high_max;

extra_variables.all_slopes = all_slopes*ones(length(history_lengths),1);