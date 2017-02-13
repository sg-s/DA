% helper function, called by tau_gain_nat_stim


function [] = plot_tau_gain_nat_stim(data,ax,c)

if nargin < 3
	c = [0 0 0];
end

deviations = data.R(:) - data.P(:);
deviations(data.R(:)<30) = NaN;
rm_this = isnan(deviations);
deviations(rm_this) = [];

all_history_lengths = [1 unique(round(logspace(1,4.5,30)))];
rho = NaN*all_history_lengths;
S = data.S(:); S = S - min(S);
for i = 1:length(all_history_lengths)
	hl = all_history_lengths(i);
	Shat = filter(ones(hl,1),hl,S);
	Shat(rm_this) = [];
	rho(i) = spear(Shat(1:10:end),deviations(1:10:end));
end



% show predictions vs. data
if ishandle(ax(1))
	x = data.P(:);
	y = data.R(:);
	r2 = oval(rsquare(x,y));
	ss = 30;
	l = plot(ax(1),x(1:ss:end),y(1:ss:end),'.','Color',c);
	legend(l,['r^2 = ' r2],'Location','southeast')
	xlabel(ax(1),'NLN model prediction (Hz)')
	ylabel(ax(1),'ab2A firing rate (Hz)')
end

if ishandle(ax(2))
	[hy,hx] = histcounts(deviations,100);
	hy = hy/sum(hy);
	hx = hx(2:end) + mean(diff(hx));
	plot(ax(2),hx,hy,'Color',c)
	xlabel(ax(2),'Deviation from NLN model (Hz)')
	ylabel(ax(2),'Probability')
end

[~,loc] = min(rho);
hl = all_history_lengths(loc);
Shat = filter(ones(hl,1),hl,S);
Shat(rm_this) = [];
[~,plot_data]=plotPieceWiseLinear(Shat,deviations,'nbins',60,'make_plot',false);

if ishandle(ax(3))
	plot(ax(3),plot_data.x,plot_data.y,'+','Color',c)
	set(ax(3),'XScale','log')
	xlabel(ax(3),['\mu_{S}  in preceding ' oval(hl) 'ms'])
	ylabel(ax(3),'Deviations from NLN model (Hz)')
end

if ishandle(ax(4))
	plot(ax(4),all_history_lengths,rho,'+-','Color',c)
	set(ax(4),'XScale','log','XTick',[1 10 100 1e3 1e4 1e5])
	xlabel(ax(4),'\tau_{gain} (ms)')
	ylabel(ax(4),'\rho')

	% also plot the stimulus and response auto-correlations
	tau_stim = autoCorrelationTime(data.S(:));
	plot(ax(4),[tau_stim tau_stim],[-1 1],':','Color',c)
	tau_resp = autoCorrelationTime(data.R(:));
	plot(ax(4),[tau_resp tau_resp],[-1 1],'--','Color',c)
	set(ax(4),'YLim',[-1 .1])
end

