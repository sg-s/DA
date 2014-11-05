allfiles= dir('*.mat');

ls_data = [];

AllControlParadigms = {};

for a = 1:length(allfiles)
	load(allfiles(a).name)
	disp(allfiles(a).name)
	for i = 1:length(spikes)
		for j = 1:width(spikes(i).A)
			disp([i j])
			if length(spikes(i).A(j,:)) > 2
				% has data, convert to f
				z = length(data(i).PID(j,:));
				time = (1:z)/SamplingRate;
				[f,t] = spiketimes2f(spikes(i).A(j,1:z),time,3e-3,3e-2);

				% which paradigm is it?
				this_paradigm = ControlParadigm(i).Name;
				this_paradigm = this_paradigm(strfind(this_paradigm,'bkg'):end);
				disp(this_paradigm)

				load_here =find(strcmp(this_paradigm, AllControlParadigms));
				if isempty(load_here)
					load_here = length(AllControlParadigms)+1;
					AllControlParadigms{load_here} = this_paradigm;
					ls_data(load_here).ORN = f(:);
					ls_data(load_here).time = t(:);
					ls_data(load_here).PID = interp1(time(:),data(i).PID(j,1:z),t(:));
				else
					ls_data(load_here).ORN = [ls_data(load_here).ORN f(:)];
					thisPID = interp1(time(:),data(i).PID(j,:),t(:));
					ls_data(load_here).PID = [ls_data(load_here).PID thisPID(:)];
				end

			end
		end
	end
end

