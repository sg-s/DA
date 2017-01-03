
% wrapper to visualise LFP slowdown models

function [] = LFP_slowdown_wrapper(p)

S = [ones(5e3,1); 2*ones(5e3,1); ones(3e3,1)];
S = repmat(S,1,10);

b = logspace(log10(.2),log10(2),10);

for i = 1:size(S,2)
	S(:,i) = S(:,i)*b(i);
end

figure('outerposition',[0 0 1501 500],'PaperUnits','points','PaperSize',[1501 500]); hold on
subplot(1,3,1); hold on
set(gca,'ColorOrder',parula(size(S,2)))
plot(S)

% solve the model
[R] = asNL(S,p);

subplot(1,3,2); hold on
set(gca,'ColorOrder',parula(size(S,2)))
plot(R)