% load spiketimes to data
load('/local-data/DA-paper/data.mat')
for i = 2:length(data)
	textbar(i,length(data))
	load(strcat('/Data/orn/data-for-da-paper/flickering-stim-carlotta/',data(i).original_name))
	data(i).spiketimes = squeeze(stA);

	% extract firing rates 
	t = 1e-4:1e-4:1e-4*length(PID);
	f = mean2(spiketimes2f(data(i).spiketimes,t,3e-3,1e-2));


	% crop data
	a =  find(data(i).full_data.time>data(i).time(1),1,'first');
	z =  find(data(i).full_data.time>data(i).time(end),1,'first');

	f = f(a:z);

	% add it to data
	data(i).ORN = f(:);

end