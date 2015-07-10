time = deltat*(1:length(ORN));

this_PID = [];
this_ORN = [];
data(load_here).PID = [];
data(load_here).ORN = [];

for i = 1:size(tracedone,1) % indexes over dilutions

	data(load_here).original_name = original_name;
	data(load_here).neuron = this_neuron;
	data(load_here).odor = this_odor;


	% calculate a time vector
	t = 0:3e-3:(deltat*length(ORN));
	this_PID = 0*t';
	this_ORN = 0*t';

	c = 1;
	for j = 1:size(tracedone,2) % indexes over repeats
		if tracedone(i,j)
			disp([i j])
			f =  spiketimes2f(squeeze(tSPKA(i,j,:)),time,3e-3,3e-2);
			this_PID = this_PID + interp1(time(:),squeeze(PID(i,j,:)),t(:));
			this_ORN = this_ORN + f;
			c=c+1;

		end
	end
	this_PID = this_PID/c;
	this_ORN = this_ORN/c;

	data(load_here).PID = [data(load_here).PID this_PID];
	data(load_here).ORN = [data(load_here).ORN this_ORN];


end 

