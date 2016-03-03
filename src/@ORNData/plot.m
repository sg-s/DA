%ORNDATA/plot.m
% overloaded plot function for the ORNDATA class
% plot makes various plots from data in the ORNDATA class
% 
% minimal usage:
% plot(orn_data,'plot_what'), e.g. plot(orn_data,'timeseries.firing_rate')
%
% where the first argument is a ORNData object, and the second argument is a 
% string specifying what sort of plot to make. 
% Use plot(orn_data,'help') to get more info.
% 
% created by Srinivas Gorur-Shandilya at 9:20 , 01 February 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [plot_handles] = plot(varargin)

% options and defaults
showr2 = false;
data_bin_type = 'pwlinear'; % can also be "dots"
nbins = 30;
ngroups = 5;
use_std = false;
plot_type = 'trial-wise'; % can also be "sem", or "mean"
grouping = [];
plot_here = [];
history_length = 500; % in ms, or the time step of the data (which is the same)
hl_min = 50; % in ms
hl_max = 1e4;
history_lengths = unique([logspace(log10(hl_min),log10(500),15) logspace(log10(500),log10(hl_max),15)]);
min_inst_gain_r2 = .8; % r2 values of inst. gain below this are discarded
min_inst_gain_firing = 10; % firing rates below 10Hz are excluded from the analysis 
show_NL = true;


% defensive programming
assert(length(varargin)>1,'Not enough input arguments.')
o = varargin{1};
assert(isa(o,'ORNData'),'1st argument should be an object from the ORNData class')
varargin(1) = [];

% figure out WHERE to plot
temp = varargin{1};
if isa(temp,'matlab.graphics.axis.Axes')
	% user has supplied a handlle to a axis; make sure it is valid
	if min(isvalid(temp))
		plot_here = varargin{1};

	else
	end
	varargin(1) = [];
else
end

% figure out WHAT to plot
allowed_plots = {'excGainAnalysis.','valveGainAnalysis.','whiffGainAnalysis.','dynamicIO.','help.','timeseries.','pdf.','Filter.','muSigma.','ioCurve.','LN.','weber.','laughlin.','instGainAnalysis.','muSigmaR2.','binnedGainAnalysis.'};
temp = varargin{1};
if isa(temp,'char')
	% user has supplied a string, make sure we can understand it 
	ok = false;
	for i = 1:length(allowed_plots)
		if ~isempty(strfind(varargin{1},allowed_plots{i}))
			ok = true;
		end
	end
	if ok 
		plot_what = varargin{1};
		varargin(1) = [];
	else
		disp('I dont know what you want me to plot. Allowed plot arguments are:')
		disp(sort(allowed_plots'))
		error('Invalid plot argument')
	end
else
	error('Expected a string specifying what to plot, got something else instead.')
end

% handle any additional arguments (they have to be options in name value syntax )
if ~isempty(varargin) 
	if iseven(length(varargin))
		for ii = 1:2:length(varargin)-1
			temp = varargin{ii};
        	if ischar(temp)
        		eval(strcat(temp,'=varargin{ii+1};'));
        	end
    	end
	else
    	error('Inputs need to be name value pairs')
    end
end

% respect use_this_segment and *se_these_trials properties
if isempty(o.use_these_trials)
	utt = 1:o.n_trials;
else
	utt = o.use_these_trials;
	if ~isempty(grouping)
		grouping = grouping(utt);
	end
	
end
if isempty(o.use_this_segment)
	uts = 1:length(o.stimulus);
else
	uts = o.use_this_segment;
end

if any(strfind(plot_what,'help'))
	disp('ORNData/plot')
	disp('You are viewing the help on how to tell plot what to plot.')
	disp('You can specify what you want to plot with a string as the 2nd or 3rd argument.')
	disp('The string comprises of three parts, separated by a period:')
	disp('plot_type.data_type.modifier')
	disp('Allowed plot types are:')
	disp(sort(allowed_plots'))
	disp('')
	disp('Allowed data types are almost any reasonable property in the ORNdata class.')
	disp('e.g., "firing_rate" is a property of the ORNData class.')
	disp('The third part, the modifier, specifies the X axis if ambiguous. e.g.,')
	disp('>> plot(orn_data,''gainAnalysis.firing.mu'')')

end

if strfind(plot_what,'Filter.')
	if isempty(plot_here)
		clear plot_here % so that plot_here gets the right class (axes object)
		figure('outerposition',[0 0 600 600],'PaperUnits','points','PaperSize',[1000 600]); hold on 			
		plot_here = gca;
	end
	% only plot fitlers. figure out which filters to plot
	if strfind(plot_what,'firing_rate')
		filtertime = o.filtertime_firing;
		K = o.K_firing(:,utt);
	elseif strfind(plot_what,'LFP')
		filtertime = o.filtertime_LFP;
		K = o.K_LFP(:,utt);
	else
		error('What to plot not specified. You told me to plot the filter, but the only things I can plot the filter of are "firing_rate" or "LFP"')
	end

	clear plot_options
	plot_options.plot_type = plot_type;
	plot_handles= plotFilters(plot_here,filtertime,K,grouping,plot_options);
	xlabel(plot_here,'Filter Lag (s)')
	ylabel(plot_here,'Filter')

end

if strfind(plot_what,'whiffGainAnalysis.')
	% verify that we get three arguments
	assert(length(strfind(plot_what,'.'))==2,'string specfiying what to plot should have 3 parts, since you are asking for a inst. gain analysis. e.g, "instGainAnalysis.firing_rate.mu')

	if isempty(plot_here)
	 	clear plot_here % so that plot_here gets the right class (axes object)
		figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on 			
		plot_here(1) = subplot(1,2,1); hold on
		plot_here(2) = subplot(1,2,2); hold on
	end

	stim = nanmean(o.stimulus(uts,utt),2);
	shat = computeSmoothedStimulus(stim,history_length);

	% figure out what to plot
	if strfind(plot_what,'firing_rate')
		pred = nanmean(o.firing_projected(uts,utt),2);
		resp = nanmean(o.firing_rate(uts,utt),2);
		ylabel(plot_here(1),'ORN Gain (Hz/V)')
	elseif strfind(plot_what,'LFP')
		pred = nanmean(o.LFP_projected(uts,utt),2);
		resp = nanmean(o.LFP(uts,utt),2);
		ylabel(plot_here(1),'ORN Gain (mV/V)')
	else
		error('What to plot not specified. You told me to plot the whiff gain analysis, but the only things I can plot the whiff gain analysis of are "firing_rate" or "LFP"')
	end
	stim = nanmean(o.firing_rate(uts,utt),2);

	% find whiffs
	[whiff_ons,whiff_offs] = findWhiffs(stim);

	% find the gains in all windows and also grab the data to plot
	[gain,gain_err,plot_data] = findGainInWindows(whiff_ons,whiff_offs,pred,resp);


	% now compute the smoothed stimulus for all history lengths and also compute the rho
	rho = NaN*history_lengths;
	for i = 1:length(history_lengths)
		temp = computeSmoothedStimulus(stim,round(history_lengths(i)));
		shat = NaN*whiff_ons;
		for j = 1:length(whiff_ons)
			shat(j) = mean(temp(whiff_ons(j):whiff_offs(j)));
		end
		rho(i) = spear(shat(gain>0 & gain_err > .8),gain(gain>0 & gain_err > .8));
	end
	plot(plot_here(2),history_lengths,rho,'k+')
	set(plot_here(2),'XScale','log')

	% show the example history length
	temp = computeSmoothedStimulus(stim,round(history_length));
	shat = NaN*whiff_ons;
	for j = 1:length(whiff_ons)
		shat(j) = mean(temp(whiff_ons(j):whiff_offs(j)));
	end
	shat = shat - nanmean(shat);
	shat = shat/nanstd(shat);
	plot(plot_here(1),shat(gain>0 & gain_err > .8),gain(gain>0 & gain_err > .8),'k+')

	% colour the plot 
	shat = shat - nanmin(shat);
	shat = shat/nanmax(shat);
	c = parula(100);
	for i = 1:length(valve_ons)
		try
			plot_handles(1).dots(i).Color = c(1+floor(shat(i)*99),:);
		catch
		end
	end
	

	% cosmetics
	xlabel(plot_here(3),'History Length (ms)')
	ylabel(plot_here(3),'\rho')
	xlabel(plot_here(2),'\mu_{stimulus} in preceding window (norm)')
	ylabel(plot_here(2),'Gain (norm)')
	xlabel(plot_here(1),'LN Model prediction (Hz)')
	ylabel(plot_here(1),'ORN Response (Hz)')
	set(plot_here(3),'YLim',[-1 1])


end

if strfind(plot_what,'excGainAnalysis.')
	% verify that we get three arguments
	assert(length(strfind(plot_what,'.'))==2,'string specfiying what to plot should have 3 parts, since you are asking for a inst. gain analysis. e.g, "instGainAnalysis.firing_rate.mu')

	if isempty(plot_here)
	 	clear plot_here % so that plot_here gets the right class (axes object)
		figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on 			
		plot_here(1) = subplot(1,3,1); hold on
		plot_here(2) = subplot(1,3,2); hold on
		plot_here(3) = subplot(1,3,3); hold on
	end

	stim = nanmean(o.stimulus(uts,utt),2);
	shat = computeSmoothedStimulus(stim,history_length);

	% figure out what to plot
	if strfind(plot_what,'firing_rate')
		pred = nanmean(o.firing_projected(uts,utt),2);
		resp = nanmean(o.firing_rate(uts,utt),2);
		ylabel(plot_here(1),'ORN Gain (Hz/V)')

		% throw out some junk points
		rm_this = pred > 2*max(resp);
		pred(rm_this) = []; resp(rm_this) = [];

	elseif strfind(plot_what,'LFP')
		pred = nanmean(o.LFP_projected(uts,utt),2);
		resp = nanmean(o.LFP(uts,utt),2);
		ylabel(plot_here(1),'ORN Gain (mV/V)')
	else
		error('What to plot not specified. You told me to plot the whiff gain analysis, but the only things I can plot the whiff gain analysis of are "firing_rate" or "LFP"')
	end

	% find all the excursions in the response
	temp = resp;
	if min(temp) >= 0
		% firing rate
	else
		% LFP
		temp = -temp;
	end
	temp = temp-min(temp);
	temp = temp/max(temp);
	[ons,offs] = computeOnsOffs(temp>.1);

	% excursions should be a certain length
	rm_this = (offs-ons) < 50 | (offs-ons) > 500;
	ons(rm_this) = []; offs(rm_this) = [];

	% find the gains in all windows and also grab the data to plot
	[gain,gain_err,plot_data] = findGainInWindows(ons,offs,pred,resp);

	if showNL
		% plot this
		for i = 1:length(plot_data)
			if gain_err(i) > .8 && ~isnan(gain_err(i)) && gain(i) > 0
				plot_handles(1).dots(i) = plot(plot_here(1),plot_data(i).x,plot_data(i).y,'k.');
			end
		end
	end
	
	% now compute the smoothed stimulus for all history lengths and also compute the rho
	rho = NaN*history_lengths;
	for i = 1:length(history_lengths)
		temp = computeSmoothedStimulus(stim,round(history_lengths(i)));
		shat = NaN*ons;
		for j = 1:length(ons)
			shat(j) = mean(temp(ons(j):offs(j)));
		end
		rho(i) = spear(shat(gain>0 & gain_err > .8),gain(gain>0 & gain_err > .8));

	end
	plot(plot_here(3),history_lengths,rho,'k+')
	set(plot_here(3),'XScale','log')

	% show the example history length
	temp = computeSmoothedStimulus(stim,round(history_length));
	shat = NaN*ons;
	for j = 1:length(ons)
		shat(j) = mean(temp(ons(j):offs(j)));
	end
	shat = shat - nanmean(shat);
	shat = shat/nanstd(shat);
	plot(plot_here(2),shat(gain>0 & gain_err > .8),gain(gain>0 & gain_err > .8),'k+')

	% colour the plot 
	shat = shat - nanmin(shat);
	shat = shat/nanmax(shat);
	c = parula(100);
	for i = 1:length(ons)
		try
			plot_handles(1).dots(i).Color = c(1+floor(shat(i)*99),:);
		catch
		end
	end
	

	% cosmetics
	xlabel(plot_here(3),'History Length (ms)')
	ylabel(plot_here(3),'\rho')
	xlabel(plot_here(2),'\mu_{stimulus} in preceding window (norm)')
	ylabel(plot_here(2),'Gain (norm)')
	xlabel(plot_here(1),'LN Model prediction (Hz)')
	ylabel(plot_here(1),'ORN Response (Hz)')
	set(plot_here(3),'YLim',[-1 1])

end


if strfind(plot_what,'valveGainAnalysis.')
	% verify that we get three arguments
	assert(length(strfind(plot_what,'.'))==2,'string specfiying what to plot should have 3 parts, since you are asking for a inst. gain analysis. e.g, "instGainAnalysis.firing_rate.mu')

	if isempty(plot_here)
	 	clear plot_here % so that plot_here gets the right class (axes object)
		figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on 			
		plot_here(1) = subplot(1,3,1); hold on
		plot_here(2) = subplot(1,3,2); hold on
		plot_here(3) = subplot(1,3,3); hold on
	end

	stim = nanmean(o.stimulus(uts,utt),2);
	shat = computeSmoothedStimulus(stim,history_length);

	% figure out what to plot
	if strfind(plot_what,'firing_rate')
		pred = nanmean(o.firing_projected(uts,utt),2);
		resp = nanmean(o.firing_rate(uts,utt),2);
		[whiff_ons,whiff_offs] = findWhiffs(pred);
		ylabel(plot_here(1),'ORN Gain (Hz/V)')
	elseif strfind(plot_what,'LFP')
		pred = nanmean(o.LFP_projected(uts,utt),2);
		resp = nanmean(o.LFP(uts,utt),2);
		temp = -pred;
		temp = temp-min(temp);
		[whiff_ons,whiff_offs] = findWhiffs(temp);
		ylabel(plot_here(1),'ORN Gain (mV/V)')
	else
		error('What to plot not specified. You told me to plot the whiff gain analysis, but the only things I can plot the whiff gain analysis of are "firing_rate" or "LFP"')
	end

	% find all valve ons and offs
	[valve_ons,valve_offs] = computeOnsOffs(o.valve(uts));
	% remove the last one because it might be at the end
	valve_ons(end) = []; valve_offs(end) = [];
	resp(resp<min_inst_gain_firing) = NaN;

	for i = 1:length(valve_ons)
		a = find(~isnan(resp(valve_ons(i):valve_offs(i))),1,'first');
		if ~isempty(a)
			a = a + valve_ons(i);
			[~,z] = max(resp(a:valve_offs(i))); z = z + a;
			if z - a > 100
				z = a + 100;
			end
			if any(isnan(resp(a:z)))
				z = a+find(~isnan(resp(a:z)),1,'last');
			end
			valve_ons(i) = a;
			try
				valve_offs(i) = z;
			catch
				valve_offs(i) = NaN;
				valve_ons(i) = NaN;
			end

		end
	end
	valve_offs(isnan(valve_offs)) = [];
	valve_ons(isnan(valve_ons)) = [];

	% find the gains in all windows and also grab the data to plot
	[gain,gain_err,plot_data] = findGainInWindows(valve_ons,valve_offs,pred,resp);

	% plot this
	if showNL
		for i = 1:length(plot_data)
			if gain_err(i) > .8 && ~isnan(gain_err(i)) && gain(i) > 0
				plot_handles(1).dots(i) = plot(plot_here(1),plot_data(i).x,plot_data(i).y,'k.');
			end
		end
	end

	% now compute the smoothed stimulus for all history lengths and also compute the rho
	rho = NaN*history_lengths;
	for i = 1:length(history_lengths)
		temp = computeSmoothedStimulus(stim,round(history_lengths(i)));
		shat = NaN*valve_ons;
		for j = 1:length(valve_ons)
			shat(j) = mean(temp(valve_ons(j):valve_offs(j)));
		end
		rho(i) = spear(shat(gain>0 & gain_err > .8),gain(gain>0 & gain_err > .8));
	end
	plot(plot_here(3),history_lengths,rho,'k+')
	set(plot_here(3),'XScale','log')

	% show the example history length
	temp = computeSmoothedStimulus(stim,round(history_length));
	shat = NaN*valve_ons;
	for j = 1:length(valve_ons)
		shat(j) = mean(temp(valve_ons(j):valve_offs(j)));
	end
	shat = shat - nanmean(shat);
	shat = shat/nanstd(shat);
	plot(plot_here(2),shat(gain>0 & gain_err > .8),gain(gain>0 & gain_err > .8),'k+')

	% colour the plot 
	shat = shat - nanmin(shat);
	shat = shat/nanmax(shat);
	c = parula(100);
	for i = 1:length(valve_ons)
		try
			plot_handles(1).dots(i).Color = c(1+floor(shat(i)*99),:);
		catch
		end
	end
	

	% cosmetics
	xlabel(plot_here(3),'History Length (ms)')
	ylabel(plot_here(3),'\rho')
	xlabel(plot_here(2),'\mu_{stimulus} in preceding window (norm)')
	ylabel(plot_here(2),'Gain (norm)')
	xlabel(plot_here(1),'LN Model prediction (Hz)')
	ylabel(plot_here(1),'ORN Response (Hz)')
	set(plot_here(3),'YLim',[-1 1])


end

if strfind(plot_what,'dynamicIO.')
 	if isempty(plot_here)
		clear plot_here % so that plot_here gets the right class (axes object)
		figure('outerposition',[0 0 600 600],'PaperUnits','points','PaperSize',[1000 600]); hold on 			
		plot_here = gca;
	end
	% only plot nonlinearity. figure out what to plot
	if strfind(plot_what,'firing_rate')
		pred = nanmean(o.firing_projected(uts,utt),2);
		stim = nanmean(o.stimulus(uts,utt),2);
		stim = computeSmoothedStimulus(stim,history_length);
		resp = nanmean(o.firing_rate(uts,utt),2);
		% respect min firing rate
		stim(resp<min_inst_gain_firing) = NaN;
		pred(resp<min_inst_gain_firing) = NaN;
		resp(resp<min_inst_gain_firing) = NaN;
		ylabel(plot_here,'Firing Rate (Hz)')

	elseif strfind(plot_what,'LFP')
		pred = nanmean(o.LFP_projected(uts,utt),2);
		resp = nanmean(o.LFP(uts,utt),2);
		ylabel(plot_here,'\DeltaLFP (mV)')
	else
		error('What to plot not specified. You told me to plot the dynamic IO curve, but the only things I can plot the dynamic ioCurve of are "firing_rate" or "LFP"')
	end
	l = labelByPercentile(nonnans(stim),ngroups);
	L = NaN*stim;
	L(~isnan(stim)) = l;
	c = parula(max(L)+1);
	if strcmp(data_bin_type,'pwlinear')
		for i = 1:max(l)
			this_pred = pred(L==i);
			this_resp = resp(L==i);
			plotPieceWiseLinear(this_pred,this_resp,'nbins',nbins,'Color',c(i,:));
		end
	elseif strcmp(data_bin_type,'dots')
		for i = 1:max(l)
			this_pred = pred(L==i);
			this_resp = resp(L==i);
			plot(plot_here,this_pred,this_resp,'.','Color',c(i,:));
		end
	end
	xlabel('Projected Stimulus (V)')
	ylabel('ORN Response (Hz)')


end

if strfind(plot_what,'ioCurve.') 
	if isempty(plot_here)
		clear plot_here % so that plot_here gets the right class (axes object)
		figure('outerposition',[0 0 600 600],'PaperUnits','points','PaperSize',[1000 600]); hold on 			
		plot_here = gca;
	end
	% only plot nonlinearity. figure out what to plot
	if strfind(plot_what,'firing_rate')
		pred = o.firing_projected(uts,utt);
		resp = o.firing_rate(uts,utt);
		ylabel(plot_here,'Firing Rate (Hz)')

	elseif strfind(plot_what,'LFP')
		pred = o.LFP_projected(uts,utt);
		resp = o.LFP(uts,utt);
		ylabel(plot_here,'\DeltaLFP (mV)')
	else
		error('What to plot not specified. You told me to plot the ioCurve, but the only things I can plot the ioCurve of are "firing_rate" or "LFP"')
	end

	clear plot_options
	plot_options.data_bin_type = data_bin_type;
	plot_options.nbins = nbins;
	plot_options.showr2 = showr2;
	plot_options.plot_type = plot_type;
	plot_options.normalise = false;
	plot_handles = plotNonlinearity(plot_here,pred,resp,grouping,plot_options);
	xlabel(plot_here,'Projected Stimulus (V)')

	if show_NL
		temp = ORNData;
		if width(resp) > 1
			resp = nanmean(resp,2);
		end
		if width(pred) > 1
			pred = nanmean(pred,2);
		end
		if strfind(plot_what,'firing_rate')
			temp.stimulus = NaN*pred;
			temp.K_firing = nanmean(o.K_firing,2);
			temp.firing_projected = pred;
			temp.firing_rate = resp;
			[~,ff] = fitNL(temp);
			l = plot(plot_here,pred,ff(pred),'r','LineWidth',2);
			% show r2
			legend(l,{['r^2 =  ' oval(rsquare(resp,ff(pred)))]})
		else
			error('not coded')
		end
	end

end

if strfind(plot_what,'LN.')
	if isempty(plot_here)
	 	clear plot_here % so that plot_here gets the right class (axes object)
		figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on 			
		plot_here(1) = subplot(1,2,1); hold on
		plot_here(2) = subplot(1,2,2); hold on
	end
	% recursively call plot with new arguments
	plot(o,plot_here(1),strrep(plot_what,'LN.','Filter.'),'plot_type',plot_type,'grouping',grouping,'data_bin_type',data_bin_type,'showr2',showr2);
	plot(o,plot_here(2),strrep(plot_what,'LN.','ioCurve.'),'plot_type',plot_type,'grouping',grouping,'data_bin_type',data_bin_type,'showr2',showr2,'nbins',nbins);

end

if strfind(plot_what,'pdf.')
	% show a distribution of ... something
	if isempty(plot_here)
	 	clear plot_here % so that plot_here gets the right class (axes object)
		figure('outerposition',[0 0 600 600],'PaperUnits','points','PaperSize',[1000 600]); hold on 			
		plot_here = gca;
	end

	if any(strfind(plot_what,'stimulus'))
		xlabel('Stimulus (V)')
		x = o.stimulus(uts,utt);
	elseif any(strfind(plot_what,'firing_rate'))
		xlabel('Firing Rate (Hz)')
		x = o.firing_rate(uts,utt);
	elseif strfind(plot_what,'firing_projected')
		xlabel('Projected Stimulus (V)')
		x = o.firing_projected(uts,utt);
	elseif any(strfind(plot_what,'LFP')) && ~any(strfind(plot_what,'LFP_projected')) && ~any(strfind(plot_what,'inst_gain'))
		xlabel('\DeltaLFP (mV)')
		x = o.LFP(uts,utt);
	elseif any(strfind(plot_what,'LFP_projected'))
		x = o.LFP_projected(uts,utt);
		xlabel('Projected Stimulus (V)')
	elseif any(strfind(plot_what,'inst_gain_firing'))
		x = o.inst_gain_firing; y = o.inst_gain_firing_err;
		r = nanmean(o.firing_rate(uts,utt),2);
		x(x<0 | y < min_inst_gain_r2 | r < min_inst_gain_firing) = [];
		x = log(x);
		xlabel(plot_here,'Inst. Gain (Hz/V)')
	elseif any(strfind(plot_what,'inst_gain_LFP'))
		x = o.inst_gain_LFP; y = o.inst_gain_LFP_err;
		x(x<0 | y < min_inst_gain_r2) = [];
		x = log(x);
		xlabel(plot_here,'Inst. Gain (mV/V)')
	elseif any(strfind(plot_what,'inst_gain_firing_err'))
		error('187 not coded')
	end

	ylabel(plot_here,'p.d.f')

	clear plot_options
	plot_options.data_bin_type = data_bin_type;
	plot_options.nbins = nbins;
	plot_options.showr2 = showr2;
	plot_options.plot_type = plot_type;
	plot_handles = plotPDF(plot_here,x,grouping,plot_options);

	if any(strfind(plot_what,'inst_gain'))
		% change the X axis to compensate for us taking the logarithm
		plot_handles.line(1).XData = exp(plot_handles.line(1).XData);
		try
			plot_handles.line(2).XData = exp(plot_handles.line(2).XData);
		catch
		end
		set(plot_here,'XSCale','log')
	end
end

if any(strfind(plot_what,'muSigma.')) && ~any(strfind(plot_what,'muSigmaR2.'))
	if isempty(plot_here)
	 	clear plot_here % so that plot_here gets the right class (axes object)
		figure('outerposition',[0 0 600 600],'PaperUnits','points','PaperSize',[1000 600]); hold on 			
		plot_here = gca;
	end

	if strfind(plot_what,'stimulus') 
		plot_this = o.stimulus(uts,utt);
		xlabel(plot_here,'\mu_{Stimulus} (V)')
		ylabel(plot_here,'\sigma_{Stimulus} (V)')
	elseif any(strfind(plot_what,'firing_rate')) && ~any(strfind(plot_what,'firing_projected'))
		plot_this = o.firing_rate(uts,utt);
		xlabel(plot_here,'\mu_{Firing Rate} (Hz)')
		ylabel(plot_here,'\sigma_{Firing Rate} (Hz)')
	elseif strfind(plot_what,'firing_projected')
		plot_this = o.firing_projected(uts,utt);
		xlabel(plot_here,'\mu_{Projected Stim.} (V)')
		ylabel(plot_here,'\sigma_{Projected Stim.} (V)')
	elseif strfind(plot_what,'LFP') && ~strfind(plot_what,'LFP_projected')
		plot_this = o.LFP(uts,utt);
	elseif strfind(plot_what,'LFP_projected')
		plot_this = o.LFP_projected(uts,utt);
	end
	plot_this = plot_this(:);

	plot_handles = plotMuSigma(plot_here,plot_this,history_length,true);
end

if any(strfind(plot_what,'muSigmaR2.'))
	if isempty(plot_here)
	 	clear plot_here % so that plot_here gets the right class (axes object)
		figure('outerposition',[0 0 600 600],'PaperUnits','points','PaperSize',[1000 600]); hold on 			
		plot_here = gca;
	end

	if strfind(plot_what,'stimulus') 
		plot_this = o.stimulus(uts,utt);
	elseif any(strfind(plot_what,'firing_rate')) && ~any(strfind(plot_what,'firing_projected'))
		plot_this = o.firing_rate(uts,utt);
	elseif strfind(plot_what,'firing_projected')
		plot_this = o.firing_projected(uts,utt);
	elseif strfind(plot_what,'LFP') && ~strfind(plot_what,'LFP_projected')
		plot_this = o.LFP(uts,utt);
	elseif strfind(plot_what,'LFP_projected')
		plot_this = o.LFP_projected(uts,utt);
	end
	ylabel(plot_here,'\rho')
	xlabel(plot_here,'History Length (ms)')
	plot_this = plot_this(:);

	r2 = NaN*history_lengths;
	for i = 1:length(history_lengths)
		[~,r2(i)] = plotMuSigma([],plot_this,history_lengths(i),false);
	end
	plot(plot_here,history_lengths,r2,'k+')
	set(plot_here,'XScale','log','YLim',[-1 1])

end

if any(strfind(plot_what,'timeseries.'))
	if isempty(plot_here)
	 	clear plot_here % so that plot_here gets the right class (axes object)
		figure('outerposition',[0 0 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold on 			
		plot_here = gca;
	end

	time = o.dt*(1:length(o.stimulus));
	if any(strfind(plot_what,'stimulus'))
		ylabel(plot_here,'Stimulus (V)')
		plot_this = o.stimulus(uts,utt);
	elseif any(strfind(plot_what,'firing_rate')) && ~any(strfind(plot_what,'firing_projected'))
		plot_this = o.firing_rate(uts,utt);
		ylabel(plot_here,'Firing Rate (Hz)')
	elseif any(strfind(plot_what,'firing_projected'))
		ylabel(plot_here,'Projected Stimulus (V)')
		plot_this = o.firing_projected(uts,utt);
	elseif any(strfind(plot_what,'LFP')) && ~any(strfind(plot_what,'LFP_projected'))
		plot_this = o.LFP(uts,utt);
		ylabel(plot_here,'\DeltaLFP (mV)')
	elseif any(strfind(plot_what,'LFP_projected'))
		plot_this = o.LFP_projected(uts,utt);
		ylabel(plot_here,'Projected Stimulus (V)')
	end
	time = time(uts);
	xlabel(plot_here,'Time (s)')

	clear plot_options
	plot_options.plot_type = plot_type;
	plot_handles = plotTimeSeries(plot_here,time,plot_this,grouping,plot_options);

end


if any(strfind(plot_what,'laughlin.'))
	if isempty(plot_here)
	 	clear plot_here % so that plot_here gets the right class (axes object)
		figure('outerposition',[0 0 400 800],'PaperUnits','points','PaperSize',[400 800]); hold on 		
		plot_here(1) = subplot(2,1,1); hold on
		plot_here(2) = subplot(2,1,2); hold on
	end
	if strfind(plot_what,'firing_rate') 
		x = o.firing_projected(uts,utt);
		y = o.firing_rate(uts,utt);
	elseif strfind(plot_what,'LFP') 
		x = o.LFP_projected(uts,utt);
		y = o.LFP(uts,utt);
	end

	% first show the PDF of the projected stimulus
	clear plot_options
	plot_options.nbins = nbins;
	plot_options.plot_type = plot_type;
	plot_handles(1).pdf = plotPDF(plot_here(1),x,grouping,plot_options);
	plot_options.data_bin_type = data_bin_type;
	plot_handles(2).h = plotLaughlin(plot_here(2),x,y,grouping,plot_options);

	% labels
	ylabel(plot_here(1),'p.d.f')
	xlabel(plot_here(2),'Projected Stimulus (V)')
	ylabel(plot_here(2),'Response (norm)')

end

if any(strfind(plot_what,'instGainAnalysis.'))

	% verify that we get three arguments
	assert(length(strfind(plot_what,'.'))==2,'string specfiying what to plot should have 3 parts, since you are asking for a inst. gain analysis. e.g, "instGainAnalysis.firing_rate.mu')

	if isempty(plot_here)
	 	clear plot_here % so that plot_here gets the right class (axes object)
		figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on 			
		plot_here(1) = subplot(1,2,1); hold on
		plot_here(2) = subplot(1,2,2); hold on
	end

	if strfind(plot_what,'.firing_rate') 
		y = o.inst_gain_firing;
		ye = o.inst_gain_firing_err;
		r = nanmean(o.firing_rate(uts,utt),2);
	elseif strfind(plot_what,'.LFP') 
		y = o.inst_gain_LFP;
		ye = o.inst_gain_LFP_err;
		r = nanmean(o.LFP(uts,utt),2);
	end
	s = nanmean(o.stimulus(uts,utt),2);

	use_mean = false;
	if strfind(plot_what,'.mu') 
		xlabel(plot_here(1),'\mu_{Stimulus} in preceding window (V)')
		use_mean = true;
	else
		xlabel(plot_here(1),'\sigma_{Stimulus} in preceding window (V)')
	end

	clear plot_options
	plot_options.use_mean = use_mean;
	plot_options.make_plot = true;
	plot_options.min_inst_gain_r2 = min_inst_gain_r2;
	plot_options.min_inst_gain_firing = min_inst_gain_firing;
	plot_options.nbins = nbins;
	plot_options.data_bin_type = data_bin_type;
	if length(history_length) > 1
		c = parula(floor(1.1*length(history_lengths))); % note that this colour map is on history_lengthS
	else
		c = [0 0 0];
	end

	for i = 1:length(history_length)
		plot_options.history_length = history_length(i);
		% find the closest hit in history_lengths
		if length(history_length) > 1
			[~,temp] = min(abs(history_lengths - history_length(i)));
		else
			temp = 1;
		end
		plot_options.colour = c(temp,:);
		plot_options.nbins = nbins;
		plot_options.use_std = use_std;
		plot_handles(1).lines(i) = plotInstGainVsStim(plot_here(1),y,ye,s,r,plot_options);
	end
	ylabel(plot_here(1),'Inst. Gain (Hz/V)')


	rho = NaN*history_lengths;
	plot_options.make_plot = false;
	plot_options.use_std = use_std;
	for i = 1:length(rho)
		plot_options.history_length = history_lengths(i);
		[~,rho(i)] = plotInstGainVsStim(plot_here(1),y,ye,s,r,plot_options);
	end

	% first plot all the history lengths in black/flat colours
	if length(history_length) == 1
		plot_handles(end+1).f2 = plot(plot_here(2),history_lengths,rho,'k+');
	else
		for i = 1:length(history_lengths)
			plot_handles(2).dots(i) = plot(plot_here(2),history_lengths(i),rho(i),'Marker','+','LineStyle','none','Color',c(i,:));
		end
	end
	set(plot_here(2),'XScale','log','YLim',[-1 1])
	xlabel(plot_here(2),'History Lengths (ms)')
	ylabel(plot_here(2),'\rho')
end


if any(strfind(plot_what,'binnedGainAnalysis.'))

	% verify that we get three arguments
	assert(length(strfind(plot_what,'.'))==2,'string specfiying what to plot should have 3 parts, since you are asking for a inst. gain analysis. e.g, "instGainAnalysis.firing_rate.mu')

	if isempty(plot_here)
	 	clear plot_here % so that plot_here gets the right class (axes object)
		figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on 			
		plot_here(1) = subplot(1,2,1); hold on
		plot_here(2) = subplot(1,2,2); hold on
	end

	if strfind(plot_what,'.firing_rate') 
		r = nanmean(o.firing_rate(uts,utt),2);
		p = nanmean(o.firing_projected(uts,utt),2);
	elseif strfind(plot_what,'.LFP') 
		r = nanmean(o.LFP(uts,utt),2);
		p = nanmean(o.LFP_projected(uts,utt),2);
	end
	s = nanmean(o.stimulus(uts,utt),2);

	use_mean = false;
	if strfind(plot_what,'.mu') 
		xlabel(plot_here(1),'\mu_{Stimulus} in preceding window (V)')
		use_mean = true;
	else
		xlabel(plot_here(1),'\sigma_{Stimulus} in preceding window (V)')
	end

	clear plot_options
	plot_options.use_mean = use_mean;
	plot_options.make_plot = true;
	plot_options.min_inst_gain_r2 = min_inst_gain_r2;
	plot_options.min_inst_gain_firing = min_inst_gain_firing;
	plot_options.nbins = nbins;
	plot_options.data_bin_type = data_bin_type;
	if length(history_length) > 1
		c = parula(floor(1.1*length(history_lengths))); % note that this colour map is on history_lengthS
	else
		c = [0 0 0];
	end

	for i = 1:length(history_length)
		plot_options.history_length = history_length(i);
		% find the closest hit in history_lengths
		if length(history_length) > 1
			[~,temp] = min(abs(history_lengths - history_length(i)));
		else
			temp = 1;
		end
		plot_options.colour = c(temp,:);
		plot_options.nbins = nbins;
		plot_options.use_std = use_std;
		plot_handles(1).lines(i) = plotBinnedGainAnalysis(plot_here(1),s,p,r,plot_options);
	end
	ylabel(plot_here(1),'Inst. Gain (Hz/V)')


	rho = NaN*history_lengths;
	plot_options.make_plot = false;
	plot_options.use_std = use_std;
	for i = 1:length(rho)
		plot_options.history_length = history_lengths(i);
		[~,rho(i)] = plotBinnedGainAnalysis(plot_here(1),s,p,r,plot_options);
	end

	% first plot all the history lengths in black/flat colours
	if length(history_length) == 1
		plot_handles(end+1).f2 = plot(plot_here(2),history_lengths,rho,'k+');
	else
		for i = 1:length(history_lengths)
			plot_handles(2).dots(i) = plot(plot_here(2),history_lengths(i),rho(i),'Marker','+','LineStyle','none','Color',c(i,:));
		end
	end
	set(plot_here(2),'XScale','log','YLim',[-1 1])
	xlabel(plot_here(2),'History Lengths (ms)')
	ylabel(plot_here(2),'\rho')
end



