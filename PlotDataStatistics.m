function [act] = PlotDataStatistics(data,td)

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
p = data(td).PID;
p = p-min(p); p =p/max(p);
[y,x] = hist(p,40);
y  = y/max(y);
plot(x,y,'k')
p = data(td).ORN;
p = p-min(p); p =p/max(p);
[y1,x1] = hist(p,40);
y1  = y1/max(y1);
plot(x1,y1,'b')
set(gca,'box','on','YLim',[0 1.5])
title('Stimulus and Response Dist.')
legend Stimulus ORN


act = [];
stop_here = 100;
while isempty(act)
	stop_here = stop_here*2;
	subplot(1,2,2), hold on
	[y,x] = autocorr(data(td).PID,stop_here);
	x=x*mean(diff(data(td).time));
	% find autocorrelation time of PID
	act = x(find(y<0,1,'first'));

end

plot(x,y,'k'), hold on
[y,x] = autocorr(data(td).ORN,stop_here);
x=x*mean(diff(data(td).time));
plot(x,y,'b'), hold on
title('Autocorrelation function of stimulus and response')
set(gca,'box','on','XLim',[min(x) max(x)])
xlabel('Lag (s)')
legend Stimulus ORN
PrettyFig;