allfiles= dir('*.mat');

dr_data = [];
c = 1;

for a = 1:length(allfiles)
	load(allfiles(a).name)
	disp(allfiles(a).name)
	for i = 1:length(spikes)
		for j = 1:width(spikes(i).A)
			disp([i j])
			if length(spikes(i).A(j,:)) > 2
				% has data, convert to f
				zz = length(data(i).PID(j,:));
				time = (1:zz)/SamplingRate;


				[f,t] = spiketimes2f(spikes(i).A(j,1:zz),time,3e-3,3e-2);


				% find out when the valve opens
				[ons,offs]=ComputeOnsOffs(ControlParadigm(i).Outputs(5,:));
				ons=floor((ons./SamplingRate)./3e-3);
				offs=ceil((offs./SamplingRate)./3e-3);
				before = 500;
				after = 800;

				PID = interp1(time(:),data(i).PID(j,:),t);
				thisPID = zeros(mean(offs-ons) + before + after+1,1);
				thisORN = thisPID;

				for k = 1:length(ons)
					thisPID = thisPID +  PID(ons(k)-before:offs(k)+after)';
					thisORN = thisORN +  f(ons(k)-before:offs(k)+after);
				end

				t = (1:length(thisPID))*3e-3;


				
				dr_data(c).ORN = thisORN;
				dr_data(c).time = t;
				dr_data(c).PID = thisPID;
				dr_data(c).divideby = length(ons);

				dr_data(c).ControlParadigmName = ControlParadigm(i).Name;

				c = c+1;
				
			end
		end
	end
end


cluster_this = [];
figure, hold on
for i =1:length(dr_data)
	plot(dr_data(i).PID)
	cluster_this = [cluster_this mean(dr_data(i).PID(end-100:end))];
	cluster_this(end) = cluster_this(end)/dr_data(i).divideby;
end

return

% cluster them
nclusters = 3; % tune the number of clusters here
id = kmeans(cluster_this(:),nclusters); 

% plot to make sure
c = parula(nclusters);
figure, hold on
for i =1:length(dr_data)
	plot(dr_data(i).time,dr_data(i).PID,'Color',c(id(i),:))
end

return
% repack the data based on this clustering



clear dr_data
dr_data = data;
save('dr_data.mat','dr_data')
