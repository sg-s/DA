%% plotMSGKinetics.m
% makes a plot LFP and firing rate delays as a function of the mean stimulus from the mean-shifted gaussians experiment

function plot_handles = plotMSGKinetics(cdata,ax)

% unpack data
v2struct(cdata)

mean_stim = mean(PID(a:z,:));

lag_LFP = NaN(1,width(PID));
lag_fA = NaN(1,width(PID));

for i = 1:width(PID)
	% compute LFP lags
	s = PID(a:z,i)-mean(PID(a:z,i)); s = s/std(s);
	x = LFP(a:z,i)-mean(LFP(a:z,i)); x = x/std(x); x = -x;
	temp = xcorr(x,s);  temp = temp/(z-a);
	[~,lag_LFP(i)] = max(temp);

	% compute firing lags
	try
		r = fA(a:z,i)-mean(fA(a:z,i)); r = r/std(r);
		temp = xcorr(r,s); temp = temp/(z-a);
		[~,lag_fA(i)] = max(temp);
	catch
	end
	
end

try
	lag_fA = lag_fA - (z-a);
	lag_fA(lag_fA<0) = NaN; lag_fA(lag_fA>1e3) = NaN;
catch
end
lag_LFP = lag_LFP - (z-a);
lag_LFP(lag_LFP<0) = NaN; lag_LFP(lag_LFP>1e3) = NaN;

c = lines(2);
[~,idx] = sort(mean_stim);
plot_handles(1) = plot(ax,mean_stim(idx),lag_LFP(idx),'+','Color',c(2,:));
xlabel(ax,'Mean Stimulus (V)')
ylabel(ax,'Lag (ms)')

plot_handles(2) = plot(ax,mean_stim(idx),lag_fA(idx),'+','Color',c(1,:));
xlabel(ax,'Mean Stimulus (V)')
set(ax,'YLim',[0 max([lag_LFP lag_fA])],'XScale','linear')

drawnow;
