%% A plot of some sine waves with different frequencies
% The plot below shows some sine waves. They have different frequencies. 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(sin(1:0.1:20),'r')
plot(sin(1:0.5:100),'k')
xlabel('Time')
ylabel('Amplitude')
legend({'Plot1','Plot2'})
title('What is this figure all about?')

%%
% You can add discussions, results, and anything else you want to say about in this homework assignment here. 

%% 
% For more help on how to use publish(), type "help publish" in MATLAB. 
