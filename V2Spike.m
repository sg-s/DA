% accepts a vector of voltages, and returns spike times for A and B
% 
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function [A,B] = V2Spike(V,dt,t_on)

% filter V 
nu=1/(2*dt); % sampling frequency 
[b,a] = butter(2,70/nu,'high');
Hd = dfilt.df2t(b,a);     % Direct-form II transposed
Vf = filter(Hd,V);         %   structure



% find peaks
p = localpeaks(Vf,'troughs');
Vp = Vf(p);
p(Vp>-3*std(Vf(1:floor(t_on/dt)))) = [];
Vp(Vp>-3*std(Vf(1:floor(t_on/dt)))) = [];


% calculate time to closest spike
td = diff(p);

% remove spurious spikes
p(td<30) = [];
td = diff(p);
td = [td(1) ; td];



% calcualte fractional amplitudes
amplitudes = repmat(Vf(p),1,20);
s = [-9:-1 1:10];
for i = 2:20
	amplitudes(:,i) = circshift(amplitudes(:,i),s(i-1));
end
clear i
frac_amp=amplitudes(:,1)./min(amplitudes')';

X = [frac_amp Vf(p)];
idx=kmeans(X,2);

plot(Vf), hold on
scatter(p(idx==2),Vf(p(idx==2)),'g')
scatter(p(idx==1),Vf(p(idx==1)),'r')

return

% extract voltage snippets
V_snip = NaN(100,length(p));
for i = 1:length(p)
	V_snip(:,i) = Vf(p(i)-50:p(i)+49);
end
clear i

% do pca
[coeff,scores]=pca(V_snip');

% cluster
X = scores(:,1:5);
idx=kmeans(X,2);

keyboard

% cluster the height of these peaks into three clusters
nh = 60;
[y,x]=hist(Vp,nh);
ff = fit(x',y','gauss3');
% discretise and sample all these distributions  
g1 = ff.a1*exp(-((x-ff.b1)/ff.c1).^2);
g2 = ff.a2*exp(-((x-ff.b2)/ff.c2).^2);
g3 = ff.a3*exp(-((x-ff.b3)/ff.c3).^2);
g = [g1; g2 ; g3];
[~,gmm]=max(g);

% find cutoff points
cp = [x(1) x(find(diff(gmm))+1) x(end)];

% split into 3 clusters
p1 = p; p1(Vp>cp(2))=[]; 
temp=logical((Vp > cp(2)).*(Vp<cp(3)));
p2 = p; p2=p2(temp);
p3 = p; p3(Vp<cp(3))=[];


% debug
plot(Vf,'k')
hold on
scatter(p1,Vf(p1),'r')
scatter(p2,Vf(p2),'b')
scatter(p3,Vf(p3),'g')

