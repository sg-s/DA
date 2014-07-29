function [act] = PlotDataStatistics(data,td,plothere)

if nargin == 2
	figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
	subplot(1,2,1), hold on
end
p = data(td).PID;
if ~isvector(p)
	% p = mean2(p);
	p = p(:);
end
p = p-min(p); p =p/max(p);
[y,x] = hist(p,50);
y  = y/max(y);
if nargin == 3
	plot(plothere(1),x,y,'k')
else
	plot(x,y,'k')
end

p = data(td).ORN;
if ~isvector(p)
	% p = mean2(p);
	p = p(:);
end
p = p-min(p); p =p/max(p);
[y1,x1] = hist(p,50);
y1  = y1/max(y1);
if nargin == 3
	plot(plothere(1),x1,y1,'b');
	set(plothere(1),'box','on','YLim',[0 1.5]);
	title(plothere(1),'Stimulus and Response Dist.')
	legend(plothere(1),{'Stimulus','ORN'});
	xlabel(plothere(1),'Normalised Value')
	ylabel(plothere(1),'Normalised Density')
else
	plot(x1,y1,'b')
	set(gca,'box','on','YLim',[0 1.5])
	title('Stimulus and Response Dist.')
	legend Stimulus ORN
end



act = [];
stop_here = 100;
p = data(td).PID;
if ~isvector(p)
	% p = mean2(p);
	p = p(:);
end
while isempty(act)
	stop_here = stop_here*2;
	[y,x] = autocorr(p,stop_here);
	x=x*mean(diff(data(td).time));
	% find autocorrelation time of PID
	act = x(find(y<0,1,'first'));

end
p = data(td).ORN;
if ~isvector(p)
	% p = mean2(p);
	p =p(:);
end
[y1,x1] = autocorr(p,stop_here);
x1=x1*mean(diff(data(td).time));

if isfield(data,'Valve')
	[y2,x2] = autocorr(data(td).Valve,stop_here);
	x2=x2*mean(diff(data(td).time));
end

if nargin == 2

	subplot(1,2,2); hold on
	plot(x,y,'k'), hold on
	plot(x1,y1,'b'), hold on
	if isfield(data,'Valve')
		plot(x2,y2,'r'), hold on
	end
	title('Autocorrelation')
	set(gca,'box','on','XLim',[min(x) max(x)])
	xlabel('Lag (s)')
	if isfield(data,'Valve')
		legend Stimulus ORN Valve
	else
		legend Stimulus ORN
	end

	PrettyFig;
else
	plot(plothere(2),x,y,'k');
	plot(plothere(2),x1,y1,'b');
	set(plothere(2),'box','on','XLim',[min(x) max(x)])
	title(plothere(2),'autocorrelation.')
	legend(plothere(2),{'Stimulus','ORN'});
	xlabel(plothere(2),'Lag (s)')
end

