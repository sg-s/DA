allfiles = dir('*.mat');
for i = 1:length(allfiles)
	load(allfiles(i).name)
	try
		spikes(8) = spikes(3);
		spikes(3).A = [];
		spikes(3).B = [];
		spikes(3).artifacts = [];
		spikes(3).amplitudes_A = [];
		spikes(3).amplitudes_B = [];
	end
	try
		spikes(3).N = [];
	end
	data(8) = data(3);
	data(3).MFC500 = [];
	data(3).PID = [];
	data(3).voltage = [];

	try
		save(allfiles(i).name,'spikes','data','-append')
	catch
		save(allfiles(i).name,'data','-append')
	end

	clearvars -except i allfiles
end
