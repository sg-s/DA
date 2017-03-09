sh = logspace(-2,1,10);
S = zeros(1e4,10);

for i = 1:size(S,2)
	S(1e3:9e3,i) = sh(i);
	S(4e3:5e3,i) = 2*sh(i);
end
S = S + 1e-3;

c = parula(11);
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
for i = 1:10
	plot(S(:,i),'Color',c(i,:))
end
set(gca,'YScale','log')
prettyFig();

