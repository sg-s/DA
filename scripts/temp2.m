example_orn = 3;

a = 1e4;
z = 5e4;

S = DNSdata.PID(a:z,DNSdata.orn==example_orn);
X = DNSdata.LFP(a:z,DNSdata.orn==example_orn);
R = DNSdata.fA(a:z,DNSdata.orn==example_orn);

[s,x,r,whiff_s,group] = binnedWhiffGainAnalysis(S,X,R,'whiff_height_frac',1/10);

figure('outerposition',[0 0 1400 901],'PaperUnits','points','PaperSize',[1400 901]); hold on
subplot(2,3,1); hold on
for i = 1:max(group)
	tempx = whiff_s(group==i)/mean(whiff_s(group==i));
	tempy = x(group==i);
	[tempx,idx] = sort(tempx);
	plot(tempx,tempy(idx),'+-');
end
l = plot(NaN,NaN,'k+');
legend(l,['\rho = ' oval(spear(whiff_s,x))],'Location','northeast')
xlabel('Whiff amplitude (norm)')
ylabel('LFP response (norm)')

subplot(2,3,2); hold on
for i = 1:max(group)
	tempx = s(group==i);
	tempy = x(group==i);
	[tempx,idx] = sort(tempx);
	plot(tempx,tempy(idx),'+-');
end
l = plot(NaN,NaN,'k+');
legend(l,['\rho = ' oval(spear(s,x))],'Location','northeast')
xlabel('\mu_{Stimulus} in preceding 300ms')
ylabel('LFP response (norm)')

subplot(2,3,4); hold on
for i = 1:max(group)
	tempx = whiff_s(group==i)/mean(whiff_s(group==i));
	tempy = r(group==i);
	[tempx,idx] = sort(tempx);
	plot(tempx,tempy(idx),'+-');
end
l = plot(NaN,NaN,'k+');
legend(l,['\rho = ' oval(spear(whiff_s,r))],'Location','northeast')
xlabel('Whiff amplitude (norm)')
ylabel('Firing response (norm)')

subplot(2,3,5); hold on
for i = 1:max(group)
	tempx = s(group==i);
	tempy = r(group==i);
	[tempx,idx] = sort(tempx);
	plot(tempx,tempy(idx),'+-');
end
l = plot(NaN,NaN,'k+');
legend(l,['\rho = ' oval(spear(s,r))],'Location','northeast')
xlabel('\mu_{Stimulus} in preceding 300ms')
ylabel('Firing response (norm)')