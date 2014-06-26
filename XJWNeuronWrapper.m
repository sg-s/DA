function [V,f,Ca] = XJWNeuronWrapper(time,stimulus,p)

stimulus = stimulus(:);
time = time(:);

nrep = 10;
noise = 0.1;


st = zeros(length(time),nrep);


for i = 1:nrep
	[V,Ca,st(:,i)] = XJWNeuronEuler(time,stimulus+noise*randn(length(time),1),p);
end

f=spiketimes2f(st,time,0.03,0.03);


if ~nargout
	figure, hold on
	plot(f)
else
	% interpolate to get back to the original length
	f = mean(f,2);
	t=time(1:10:end);
	f = interp1(t,f,time);
end 