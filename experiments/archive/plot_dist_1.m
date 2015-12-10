figure('outerposition',[0 0 1000 800],'PaperUnits','points','PaperSize',[1000 800]); hold on
subplot(2,2,1), hold on
c = jet(length(data));
use_these = 2:length(data)-1;
for i = use_these
        for j = 1:width(data(i).PID)
                [y,x] = hist(data(i).MFC500(j,200000:500000),50);
                plot(x,y,'Color',c(i,:))
        end
end
xlabel('MFC Flow (V)')
ylabel('Count')
subplot(2,2,2), hold on
for i = use_these
        for j = 1:width(data(i).PID)
                [y,x] = hist(data(i).PID(j,200000:500000),50);
                plot(x,y,'Color',c(i,:))
        end
end
xlabel('PID (V)')
ylabel('Count')
subplot(2,2,3:4), hold on
for i = use_these
        if ~isempty(data(i).PID)
                time = 1e-4*(1:length(data(i).PID));
                plot(time,mean(data(i).PID),'Color',c(i,:))
        end
end
set(gca,'XLim',[20 30])
xlabel('Time (s)')
ylabel('PID (V)')