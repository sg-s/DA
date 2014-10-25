allfiles= dir('*.mat');

ls_data = [];

AllControlParadigms = {};

for a = 1:length(allfiles)
	load(allfiles(a).name)
	for i = 1:length(spikes)
		for j = 1:width(spikes(i).A)
			disp([i j])
			if length(spikes(i).A(j,:)) > 2
				% has data, convert to f
				time = (1:length(spikes(i).A(j,:)))/SamplingRate;
				[f,t] = spiketimes2f(spikes(i).A(j,:),time,3e-3,3e-2);
				load_here =find(strcmp(ControlParadigm(i).Name, AllControlParadigms));
				if isempty(load_here)
					load_here = length(AllControlParadigms)+1;
					AllControlParadigms{load_here} = ControlParadigm(i).Name;
					ls_data(load_here).ORN = f(:);
					ls_data(load_here).time = t(:);
					ls_data(load_here).PID = interp1(time(:),data(i).PID(j,:),t(:));
				else
					ls_data(load_here).ORN = [ls_data(load_here).ORN f(:)];
					thisPID = interp1(time(:),data(i).PID(j,:),t(:));
					ls_data(load_here).PID = [ls_data(load_here).PID thisPID(:)];
				end

			end
		end
	end
end