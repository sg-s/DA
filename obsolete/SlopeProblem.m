%% Why is the Slope not 1 for Damon's function?
% The best fit line to the prediction from a linear filter and the data should be 1. But it isn't. 
% 
%%
% First we make an Gaussian noise input and create an exponential filter, and make the output by convolving the filter with the input.
filter_length = 333;
a = rand(1,10000);
a = PID;
filtertime = 1e-3:1e-3:filter_length*1e-3;
Kexp = exp(-50*filtertime);
f = filter(Kexp,1,a);

%%
% We can now reconstruct the filter:
K = Back_out_1dfilter_new(a,f,filter_length,0);
fp = filter(K,1,a);
figure, hold on
plot(Kexp), hold on
plot(K(1:end-1),'r')

%% 
% and the slope of best fit is:
[fall gof] = fit(fp(filter_length+2:end)',f(filter_length+2:end)','Poly1');
disp(fall.p1)

%%
% Now, we add some Gaussian noise to the output. 
f = f + randn(1,length(a));
K = Back_out_1dfilter_new(a,f,filter_length,0);
fp = filter(K,1,a);
figure, hold on
plot(Kexp), hold on
plot(K(1:end-1),'r')

%% 
% and the slope of best fit is:
[fall gof] = fit(fp(filter_length+2:end)',f(filter_length+2:end)','Poly1');
disp(fall.p1)

