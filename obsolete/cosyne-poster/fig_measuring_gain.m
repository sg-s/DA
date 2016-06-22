

% make figure placeholders 
fig_handle=figure('outerposition',[0 0 800 600],'PaperUnits','points','PaperSize',[800 600]); hold on
clf(fig_handle);

for i = 1:4
	axes_handles(i) = subplot(2,2,i);
	hold(axes_handles(i),'on');
end


[PID, LFP, fA, paradigm, orn, fly, AllControlParadigms, paradigm_hashes] = consolidateData('/local-data/DA-paper/LFP-MSG/september',1);

% remove baseline from all PIDs
for i = 1:width(PID)
	PID(:,i) = PID(:,i) - mean(PID(1:5e3,i));
end

% sort the paradigms sensibly
sort_value = [];
for i = 1:length(AllControlParadigms)
	sort_value(i) = (mean(AllControlParadigms(i).Outputs(1,:)));
end
[~,idx] = sort(sort_value);


AllControlParadigms = AllControlParadigms(idx);
paradigm_new = paradigm*NaN;
for i = 1:length(idx)
	paradigm_new(paradigm == idx(i)) = i;
end
paradigm = paradigm_new;

% throw our bad traces
bad_trials = (sum(fA) == 0 | isnan(sum(fA)));
PID(:,bad_trials) = [];
LFP(:,bad_trials) = [];
fA(:,bad_trials) = [];
paradigm(bad_trials) = [];
orn(bad_trials) = [];

% extract filters and find gain
a = 35e3; z = 55e3;
[K,fA_pred,fA_gain] = extractFilters(PID,fA,'use_cache',true,'a',a,'z',z);

% some core variables
dt = 1e-3;
c = parula(max(paradigm)+1); % colour scheme

% plot some stimulus
example_dose = 1;
plot_this = nanmean(PID(:,paradigm==example_dose),2);
time = dt*(1:length(plot_this));
axes(axes_handles(1))
plot(time,plot_this,'Color',c(example_dose,:));
set(axes_handles(1),'XLim',[40 50],'YLim',[0 .55])


% plot the filter for the lowest dose 
filtertime = (-100:599)*1e-3;
axes(axes_handles(2))
shadedErrorBar(filtertime,mean(K(:,paradigm==example_dose),2),sem(K'),{'Color',c(example_dose,:)})
set(axes_handles(2),'XLim',[-.1 .6])

% plot linear prediction vs. data for the example dose. 
y = mean(fA(a:z,paradigm==example_dose),2);
x = mean(fA_pred(a:z,paradigm==example_dose),2);
time = 40+dt*(1:length(x));

[ax,plot1,plot2] = plotyy(axes_handles(3),time,y,time,x);
set(ax(1),'XLim',[35 55],'YLim',[min(y(5e3:end)) max(y(5e3:end))],'YColor',c(example_dose,:))
set(ax(2),'XLim',[35 55],'YLim',[min(x(5e3:end)) max(x(5e3:end))],'YColor','r')
set(plot1,'Color',c(example_dose,:))
set(plot2,'Color','r')

set(axes_handles(2),'box','off')
set(ax(2),'XMinorTick','on','YMinorTick','on','YTick',[0:0.1:0.5],'YTickLabel',{'0','.1','.2','.3','.4','.5'},'XLim',[40 50])
set(ax(1),'XMinorTick','on','YMinorTick','on','YTick',[0:20:60],'YTickLabel',{'0','20','40','60'},'XLim',[40 50])
set(axes_handles(3),'box','off')

plot(axes_handles(4),x(1:25:end),y(1:25:end),'.','Color',c(example_dose,:));
ff= fit(x(~isnan(x)),y(~isnan(x)),'poly1');
l = plot(axes_handles(4),sort(x),ff(sort(x)),'k');
L = {};
L{1} = strcat('Gain=',oval(ff.p1),'Hz/V, r^2=',oval(rsquare(y,x)));


set(axes_handles(4),'YColor',c(example_dose,:),'XColor','r')
set(axes_handles(3:4),'YLim',[0 65])

movePlot(axes_handles(2),'right',.08)
movePlot(axes_handles(4),'right',.08)

prettyFig('plw=2;','lw=1.5;','fs=18;')



