% fit all data
load('/local-data/DA-paper/data.mat')

NLNParam = zeros(4,21);

for i = 2:21
	[~,x] = FitNLNModel(d,[],.99);
	NLNParam(:,i) = x;
end