% fit all data
load('/local-data/DA-paper/data.mat')

NLNParam = zeros(4,21);

for i = 20:21
	clear d
	d.stimulus = data(i).PID;
	d.stimulus(d.stimulus<0) = 0;
	d.response = data(i).ORN;
	[~,x] = FitNLNModel(d,[],.99);
	NLNParam(:,i) = x;
end