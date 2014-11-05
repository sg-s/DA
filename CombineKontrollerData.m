% CombineKontrollerData.m
% combines data from mulitple kontroller-generated data files (defaults to all .mats) and automatically and intelligently combines them
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

allfiles= dir('*.mat');

% define what you want your output to be called here
combined_data = [];

AllControlParadigms = {};

for a = 1:length(allfiles)
	load(allfiles(a).name)
	for i = 1:length(spikes)
		for j = 1:width(spikes(i).A)
			disp([i j])
			if length(spikes(i).A(j,:)) > 2
				% has data, convert to f
				z = length(data(i).PID(j,:));
				time = (1:z)/SamplingRate;
				[f,t] = spiketimes2f(spikes(i).A(j,1:z),time,3e-3,3e-2);

				% define the category it is in here
				pname = ControlParadigm(i).Name;
				pname = pname(strfind(pname,'_bkg'):end);
				load_here =find(strcmp(pname, AllControlParadigms));


				if isempty(load_here)
					load_here = length(AllControlParadigms)+1;
					AllControlParadigms{load_here} = pname;
					combined_data(load_here).ORN = f(:);
					combined_data(load_here).time = t(:);
					combined_data(load_here).PID = interp1(time(:),data(i).PID(j,1:z),t(:));
				else
					combined_data(load_here).ORN = [combined_data(load_here).ORN f(:)];
					thisPID = interp1(time(:),data(i).PID(j,:),t(:));
					combined_data(load_here).PID = [combined_data(load_here).PID thisPID(:)];
				end

			end
		end
	end
end