% accepts a vector of voltages, and returns spike times for A and B
% 
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function [A,B] = V2Spike(V,dt,t_on)
V = V(:);
time = dt:dt:dt*length(V);

% two spikes cannot be closer than 1.4ms

% filter V 
nu=1/(2*dt); % sampling frequency 
[b,a] = butter(2,70/nu,'high');
Hd = dfilt.df2t(b,a);     % Direct-form II transposed
Vf = filter(Hd,V);         %   structure


% find peaks
p = find(localpeaks(Vf,'troughs'));
Vp = Vf(p);


% throw out all noise peaks
p(Vp>-2*std(Vf(1:floor(t_on/dt)))) = [];
Vp = Vf(p);

% throw out outliers
p(Vp<10*mean(Vp)) = [];
Vp = Vf(p);

% calculate time to closest spike
td = diff(p);

% remove spurious spikes -- 1: resolve almost co-incident spikes
resolve_these = find(td<10);
conflicting_pairs = [resolve_these (resolve_these+1)];

remove_these = [];
for i = 1:length(conflicting_pairs)
	if Vf(p(conflicting_pairs(i,1))) > Vf(p(conflicting_pairs(i,2)))
		remove_these = [remove_these conflicting_pairs(i,1)];
	else
		remove_these = [remove_these conflicting_pairs(i,2)];
	end
	
end
clear i conflicting_pairs resolve_these
p(remove_these) = [];
td = diff(p);

% remove spurious spikes -- 2: remove shallow minima
remove_these = [];
for i = 1:length(p)
	if p(i) > 2 && p(i)+2 < length(Vf) 
		if Vf(p(i)) <  Vf(p(i)-2) && Vf(p(i)+2)
		else
			remove_these = [remove_these i];
		end
	end
	
end
clear i 
p(remove_these) = [];

% we can't do spikes in the first 30 samples
p(p<31) = [];


% find preceding maxima for each spike
max_amp = NaN*p;
max_amp_loc = NaN*p;
for i = 1:length(max_amp)
	this_snippet = Vf(p(i)-30:p(i));
	[max_amp(i) max_amp_loc(i)]=max(this_snippet);
	max_amp_loc(i) = max_amp_loc(i) + p(i) -31;
end
clear i


% % build a density estimate of all spikes
% d = 0*V;
% w = 500;
% for i = 1:length(p)
% 	try
% 		x = p(i)-w:p(i)+w;
% 		b = p(i);
% 		c = 500;
% 		d(p(i)-w:p(i)+w) = d(p(i)-w:p(i)+w)' + exp(-((x-b)/c).^2);
% 	catch
% 	end
% end



% correct amplitudes only if density estimate is above some threshold
d_thresh = 4;
amplitudes = max_amp - Vf(p);
amplitudes = repmat(amplitudes,1,20);

% calculate fractional amplitudes as a fraction of the largest amplitude of spike
s = [-9:-1 1:10];
% s = [-19:-1];
for i = 2:20
	amplitudes(:,i) = circshift(amplitudes(:,i),s(i-1));
end
clear i
frac_amp=amplitudes(:,1)./max(amplitudes(:,2:20)')';
frac_amp2=amplitudes(:,1)./min(amplitudes(:,2:20)')';


% for i = 1:length(p)
% 	% find amplitudes of all spikes in the preceding 1s
% 	prev_spikes = find((p>(p(i)-10000)).*(p<p(i)));
% 	if isempty(prev_spikes)
% 		% ok, let's use all spikes in the next 100ms instead
% 		prev_spikes = find((p<(p(i)+10000)).*(p>p(i)));
% 	end
% 	keyboard
% end


% fit two gaussians
nh = 40;
[y,x]=hist(frac_amp,nh);
ff = fit(x',y','gauss2');
% discretise and sample all these distributions  
g1 = ff.a1*exp(-((x-ff.b1)/ff.c1).^2);
g2 = ff.a2*exp(-((x-ff.b2)/ff.c2).^2);
g = [g1; g2];
[~,gmm]=max(g);

% find cutoff point between the two gaussians
cutoff = x(max((find(diff(gmm))+1)));
p1 = p(frac_amp<cutoff);
p2 = p(frac_amp>cutoff);

frac_amp1 = frac_amp(frac_amp<cutoff);
frac_amp2 = frac_amp(frac_amp>cutoff);

% assign to A and B neuron
if mean(Vf(p1)) < mean(Vf(p2))
	A = p1; B = p2;
else
	A = p2; B = p1;
end


if nargout == 0
	plot(time,Vf,'k'), hold on
	scatter(time(p1),Vf(p1),'b')
	scatter(time(p2),Vf(p2),'r')
	% for i = 1:length(p1)
	% 	text(p1(i),Vf(p1(i))-0.01,oval(frac_amp1(i),2),'Color','b')
	% end
	% for i = 1:length(p2)
	% 	text(p2(i),Vf(p2(i))-0.01,oval(frac_amp2(i),2),'Color','r')
	% end

	scatter(time(max_amp_loc),max_amp,'kx')

end

return




