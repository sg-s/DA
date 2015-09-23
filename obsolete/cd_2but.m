% temp file that loads data from carlotta's experiments
load('/local-data/DA-paper/raw-data/carlotta-fig-4/data_2011_11_02_ab3A_2but_noback_N8_SPK.mat')
% wild assumptions
T = 11;
time = 1e-4:1e-4:T;
dt = 3e-3;
time2 = 0:3e-3:T;

% get the firing rates
f = NaN(size(tSPKA,1),length(time2));

for i = 1:size(tSPKA,1)
	% for each conc
	ftemp = spiketimes2f(squeeze(tSPKA(i,:,:)),time,3e-3,3e-2);
	ftemp(:,max(ftemp)==0) = []; 
	f(i,:) = mean2(ftemp);
end



load('/local-data/DA-paper/raw-data/carlotta-fig-4/data_2011_11_02_ab3A_2but_back2_N8_SPK.mat')
% get the firing rates
f2 = NaN(size(tSPKA,1),length(time2));

for i = 1:size(tSPKA,1)
	% for each conc
	ftemp = spiketimes2f(squeeze(tSPKA(i,:,:)),time,3e-3,3e-2);
	ftemp(:,max(ftemp)==0) = []; 
	if ~isvector(ftemp)
		f2(i,:) = mean2(ftemp);
	else
		f2(i,:) = ftemp;
	end
end
f(12:13,:) = [];
f2(12:13,:) = [];

figure, hold on
subplot(1,2,1), hold on

c = jet(size(f,1));
for i = 1:size(f,1)
	plot(time2,f(i,:),'Color',c(i,:))
end

subplot(1,2,2), hold on

for i = 1:size(f,1)
	plot(time2,f2(i,:),'Color',c(i,:))
end

subplot(1,3,3), hold on
plot(max(f'),'k')
plot(max(f2'),'r')
