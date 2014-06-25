function [f] = XJWNeuronWrapper(time,stimulus,p)

nrep = 20;
noise = 0.1;
dt = mean(diff(time));

st = zeros(length(time),nrep);

for i = 1:nrep
	st(:,i) = XJWNeuronEuler(time,stimulus+noise*randn(length(time),1),p);
end

f=(mean(spiketimes2f(st,time,10*dt)));