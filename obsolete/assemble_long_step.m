allfiles= dir('*.mat');

ls_data = [];
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



				% find out when the background valve opens
				aa=find(ControlParadigm(i).Outputs(6,:),1,'first');
				aa=ceil((aa/SamplingRate)/3e-3);
				if isempty(aa)
					aa = 1;
				end

				% only register cases where we have some blank
				if aa > 501
					aa = aa - 500;

					[f,t] = spiketimes2f(spikes(i).A(j,1:zz),time,3e-3,3e-2);


					% find out when the valve opens
					[ons,offs]=ComputeOnsOffs(ControlParadigm(i).Outputs(5,:));
					z=floor((ons(1)/SamplingRate)/3e-3);
					ls_data(c).ORN = f(aa:z);
					ls_data(c).time = t(aa:z);
					ls_data(c).PID = interp1(time(:),data(i).PID(j,:),t(aa:z));

					ls_data(c).ControlParadigmName = ControlParadigm(i).Name;

					figure, hold on
					plot(ls_data(c).time,ls_data(c).PID)
					title(mat2str(c))
					pause(1)
					delete(gcf)
					
				

					c = c+1;
				end
			end
		end
	end
end

cluster_this = [];
figure, hold on
for i =1:length(ls_data)
	plot(ls_data(i).time,ls_data(i).PID)
	cluster_this = [cluster_this mean(ls_data(i).PID(end-100:end))];
end


% cluster them
nclusters = 2; % tune the number of clusters here
id = kmeans(cluster_this(:),nclusters); 

% plot to make sure
c = jet(nclusters);
figure, hold on
for i =1:length(ls_data)
	plot(ls_data(i).time,ls_data(i).PID,'Color',c(id(i),:))
end


% repack the data based on this clustering
clear data
for i = 1:nclusters
	data(i).PID = [];
	data(i).ORN = [];
	do_these = find(id==i);
	for j = 1:length(do_these)
		do_this = do_these(j);
		data(i).PID = [data(i).PID ;ls_data(do_this).PID];
		data(i).ORN = [data(i).ORN ;ls_data(do_this).ORN'];
	end
end

clear ls_data
ls_data = data;
save('ls_data.mat','ls_data')
