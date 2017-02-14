% helper function, called by tau_gain_nat_stim


function [rho,all_history_lengths] = plot_tau_gain_nat_stim(data,ax,varargin)

% options and defaults
options.c = [0 0 0];
options.response_cutoff = 30;
options.method = 'Spearman';

if nargout && ~nargin 
	varargout{1} = options;
    return
end

% validate and accept options
if iseven(length(varargin))
	for ii = 1:2:length(varargin)-1
	temp = varargin{ii};
    if ischar(temp)
    	if ~any(find(strcmp(temp,fieldnames(options))))
    		disp(['Unknown option: ' temp])
    		disp('The allowed options are:')
    		disp(fieldnames(options))
    		error('UNKNOWN OPTION')
    	else
    		options.(temp) = varargin{ii+1};
    	end
    end
end
elseif isstruct(varargin{1})
	% should be OK...
	options = varargin{1};
else
	error('Inputs need to be name value pairs')
end

deviations = data.R(:) - data.P(:);
S = data.S(:); S = S - min(S);
deviations(data.R(:)<options.response_cutoff) = NaN;

all_history_lengths = [1 unique(round(logspace(1,4.5,30)))];
rho = NaN*all_history_lengths;

for i = 1:length(all_history_lengths)
	hl = all_history_lengths(i);
	Shat = filter(ones(hl,1),hl,S);


	% remove nans
	rm_this = isnan(Shat) | isnan(deviations);
	x = Shat(~rm_this);
	y = deviations(~rm_this);

	switch options.method
	case 'Spearman'
		rho(i) = corr(x(1:10:end),y(1:10:end),'type','Spearman');
	case 'Pearson'
		rho(i) = corr(x(1:10:end),y(1:10:end),'type','Pearson');
	case 'slope'
		temp = fit(x(1:10:end),y(1:10:end),'poly1');
		rho(i) = temp.p1;
	end
end

% show predictions vs. data
if ishandle(ax(1))
	x = data.P(:);
	y = data.R(:);
	r2 = oval(rsquare(x,y));
	ss = 30;
	l = plot(ax(1),x(1:ss:end),y(1:ss:end),'.','Color',options.c);
	legend(l,['r^2 = ' r2],'Location','southeast')
	xlabel(ax(1),'NLN model prediction (Hz)')
	ylabel(ax(1),'ab2A firing rate (Hz)')
end

if ishandle(ax(2))
	[hy,hx] = histcounts(deviations,100);
	hy = hy/sum(hy);
	hx = hx(2:end) + mean(diff(hx));
	plot(ax(2),hx,hy,'Color',options.c)
	xlabel(ax(2),'Deviation from NLN model (Hz)')
	ylabel(ax(2),'Probability')
end



if ishandle(ax(3))

	[~,loc] = min(rho);
	hl = all_history_lengths(loc);
	Shat = filter(ones(hl,1),hl,S);
	rm_this = isnan(deviations) | isnan(Shat);

	[~,plot_data]=plotPieceWiseLinear(Shat(~rm_this),deviations(~rm_this),'nbins',60,'make_plot',false);

	plot(ax(3),plot_data.x,plot_data.y,'+','Color',options.c)
	set(ax(3),'XScale','log')
	xlabel(ax(3),['\mu_{S}  in preceding ' oval(hl) 'ms'])
	ylabel(ax(3),'Deviations from NLN model (Hz)')
end

if ishandle(ax(4))
	plot(ax(4),all_history_lengths,rho,'+-','Color',options.c)
	set(ax(4),'XScale','log','XTick',[1 10 100 1e3 1e4 1e5])
	xlabel(ax(4),'\tau_{gain} (ms)')
	switch options.method
	case 'Spearman'
		ylabel(ax(4),'Spearman''s \rho')
	case('Pearson')
		ylabel(ax(4),'Pearson''s r')
	case('slope')
		ylabel(ax(4),'Slope (Hz/V)')
	end

	% also plot the stimulus and response auto-correlations
	tau_stim = autoCorrelationTime(data.S(:));
	plot(ax(4),[tau_stim tau_stim],[-1 1],':','Color',options.c)
	tau_resp = autoCorrelationTime(data.R(:));
	plot(ax(4),[tau_resp tau_resp],[-1 1],'--','Color',options.c)
	set(ax(4),'YLim',[-1 .1])
end

