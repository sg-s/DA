
%  ######         ###       ####    ##    ##                   
% ##    ##       ## ##       ##     ###   ##                   
% ##            ##   ##      ##     ####  ##                   
% ##   ####    ##     ##     ##     ## ## ##                   
% ##    ##     #########     ##     ##  ####                   
% ##    ##     ##     ##     ##     ##   ###                   
%  ######      ##     ##    ####    ##    ##   


%    ###       ########        ###       ########     ######## 
%   ## ##      ##     ##      ## ##      ##     ##       ##    
%  ##   ##     ##     ##     ##   ##     ##     ##       ##    
% ##     ##    ##     ##    ##     ##    ########        ##    
% #########    ##     ##    #########    ##              ##    
% ##     ##    ##     ##    ##     ##    ##              ##    
% ##     ##    ########     ##     ##    ##              ##    

%%
% Do ORNs rapidly modulate gain?

% Now that we know that we can fit the model, and that our optimisation algorithm works, we will try to fit real data. The following figure shows the sample data we use. The top panel shows the stimulus as measured by the PID, and the bottom panel shows the firing rate of the ORN. The firing rate has been divided by its standard deviation and has been mean subtracted. Also shown is the linear prediction of the firing rates. 


if ~exist('PID','var')
	filename = '~/Desktop/final_2011_06_14_ab3A_1o3ol3X-3_20ml_30sec_30ms_rand.mat';
	[PID, time, f,Valve,uncropped] = PrepData3(filename);
	PID = PID(:);
	time = time(:);
	f = f(:);

	% detrend PID
	ptrend = fit(time,PID,'Poly1'); 
	PID = PID - (ptrend(time) - mean(ptrend(time)));

	% prepare data
	f = f/std(f);
	f = f(:) - mean(f);


	% assemble into a data structure
	data.PID = PID;
	data.f = f;
	data.Valve = Valve;
end


% build a simple linear model
K = FindBestFilter(PID(500:end),f(500:end),[],'min_cutoff=min(response);');
LinearFit = filter(K,1,PID-mean(PID))+ mean(f(500:end));
LinearFit(LinearFit<min(f)) = min(f);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
a=multiplot(time,PID,f,LinearFit);
title(a(1),'Figure 2: ORN data, and linear model prediction')
set(gca,'XLim',[20 25])
PrettyFig;


%%
% Now, we fit the DA model to the data using a pattern search optimisation (the type of optimisation doesn't matter. Pattern search is faster and converges to local minima faster than GA). $\gamma$ is constrained to $[0,1]$

% use pattern search to find parameters
if ~exist('psp.mat','file')
	x0 = [3557  203   0.35 1.8 1.9   17   1];
	lb = [1200  30    0    0   1     0    1];
	ub = [5200  300   1    10  5     30   5];
	psoptions = psoptimset('UseParallel',true,'CompletePoll', 'on', 'Vectorized', 'off','Display','iter','MaxIter',2000,'MaxFunEvals',10000);
	x = patternsearch(@(x) DA_cost_function(x,data,@Cost,'ga'),x0,[],[],[],[],lb,ub,psoptions);
	psp = ValidateDAParameters(x,'ga');
	save('psp.mat','psp')
end

load psp.mat
figure('outerposition',[0 0 1400 500],'PaperUnits','points','PaperSize',[1400 500]); hold on
DAFit = DA_integrate(PID,psp);
DAFit = DAFit - mean(DAFit(500:end));
DAFit(DAFit<min(f))=min(f);

mpo = [];
mpo.legend = 1;

subplot(1,2,1), hold on
multiplot(time,f,LinearFit,DAFit,mpo);
set(gca,'XLim',[min(time)+1 min(time)+3])
PrettyFig;

mpo.legend = 0;

subplot(1,2,2),  hold on
multiplot(time,f,LinearFit,DAFit,mpo);
set(gca,'XLim',[mean(time)-1 mean(time)+1])
set(gca,'YTick',[])
set(gca,'YColor','w')
PrettyFig;






%%
% The DA model seems to do a pretty good job estimating the ORN output. Is it better than the simple linear prediction? Here, we compare the r-square and the l-2 norm between the data and the fit. For the simple linear model, the r-square is 
s = 500;
z = length(f);
disp(rsquare(f(s:end),LinearFit(s:end)))

%%
% cf. for the DA model fit, the rsquare is 
disp(rsquare(f(s:end),DAFit(s:end)))

%% 
% The l-2 norm of the linear fit is
disp(l2(f(s:end),LinearFit(s:end)))

%% 
% cf. l-2 norm of the DA fit is
disp(l2(f(s:end),DAFit(s:end)))



%% Real Data: Gain Analysis
% We can perform a similar gain analysis like we did on the synthetic data on the real data from the ORN. The plot on the left compares the ORN response data (on the Y-axis) to the linear fit, while the plot on the right compares the data to the DA model fit. 
s = 300; % when we start for the gain analysis
z = length(f); % where we end
history_lengths = [0.102];
% hl = history_lengths/3; % history lengths better be divisible by 3!
% shat = NaN(length(hl),length(PID(s:z)));
% for i = 1:length(hl)
% 	shat(i,:) = filtfilt(ones(1,hl(i))/hl(i),1,PID(s:z));
% 	shat(i,1:hl(i)) = NaN;
% end

x.data = f(s:z);
x.prediction = LinearFit(s:z);
x.stimulus = PID(s:z);
x.time = time(s:z);

% make gain analysis plot for synthetic data and linear model
figure('outerposition',[0 0 900 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plothere=subplot(1,2,1); hold on
GainAnalysis2(x,history_lengths,filter_length,'plotid',1,'plothere',plothere);
xlabel('Linear Prediction','FontSize',font_size)
ylabel('ORN response (a.u.)','FontSize',font_size)

% make gain analysis plot for synthetic data and DA model
x.data = f(s:z);
x.prediction = DAFit(s:z);
x.stimulus = PID(s:z);
x.time = time(s:z);
plot_here=subplot(1,2,2); hold on
GainAnalysis2(x,history_lengths,filter_length,'plotid',1,'plothere',plothere);
xlabel('DA Prediction','FontSize',font_size)
legend('Location',[0.7674    0.2927    0.21    0.1370],{'all data','bottom 10%','top 10%'})



%%
% In the analysis above, we have kept the "window history" length fixed at ~100ms. How does varying this window change the separation of slopes of best fit lines for the top 10% and the bottom 10%? 

%%
% And we can do the same thing for the ORN response data.

s = 300; % when we start for the gain analysis
history_lengths = [.030 .102 .150 .300 0.600 1.002 1.500 2.001];

x.data = f(s:z);
x.prediction = LinearFit(s:z);
x.stimulus = PID(s:z);
x.time = time(s:z);


figure('outerposition',[0 0 900 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
% make gain analysis plot for synthetic data and linear model
plothere=subplot(1,2,1); hold on
GainAnalysis2(x,history_lengths,filter_length,'plotid',2,'plothere',plothere);
title('Linear Prediction','FontSize',font_size)
set(gca,'YLim',[0.7 1.5])

x.prediction = DAFit(s:z);
% make gain analysis plot for synthetic data and DA model
plothere=subplot(1,2,2); hold on
GainAnalysis2(x,history_lengths,filter_length,'plotid',2,'plothere',plothere);
ylabel('ORN response (a.u.)','FontSize',font_size)
legend('Location',[0.7674    0.5927    0.21    0.1370],{'all data','bottom 10%','top 10%'})
set(gca,'YLim',[0.7 1.5])






% ########        ###                                         
% ##     ##      ## ##                                        
% ##     ##     ##   ##                                       
% ##     ##    ##     ##                                      
% ##     ##    #########                                      
% ##     ##    ##     ##                                      
% ########     ##     ##              


% ##     ##     #######     ########     ########    ##       
% ###   ###    ##     ##    ##     ##    ##          ##       
% #### ####    ##     ##    ##     ##    ##          ##       
% ## ### ##    ##     ##    ##     ##    ######      ##       
% ##     ##    ##     ##    ##     ##    ##          ##       
% ##     ##    ##     ##    ##     ##    ##          ##       
% ##     ##     #######     ########     ########    ######## 


% ########    ####    ########                                
% ##           ##        ##                                   
% ##           ##        ##                                   
% ######       ##        ##                                   
% ##           ##        ##                                   
% ##           ##        ##                                   
% ##          ####       ##    

%%
% Can the DA Model account for fast adaptation? 












% ########     ####    ########    ########        
% ##     ##     ##     ##          ##              
% ##     ##     ##     ##          ##              
% ##     ##     ##     ######      ######          
% ##     ##     ##     ##          ##              
% ##     ##     ##     ##          ##              
% ########     ####    ##          ##   


%  #######  ########   #######  ########   ######  
% ##     ## ##     ## ##     ## ##     ## ##    ## 
% ##     ## ##     ## ##     ## ##     ## ##       
% ##     ## ##     ## ##     ## ########   ######  
% ##     ## ##     ## ##     ## ##   ##         ## 
% ##     ## ##     ## ##     ## ##    ##  ##    ## 
%  #######  ########   #######  ##     ##  ######     

%%
% How does fast adaptation depend on odor type?






%  #######  ########   #######  ##     ## ########        
% ##     ## ##     ## ##     ## ##     ## ##     ##       
% ##     ## ##     ## ##     ## ##     ## ##     ##       
% ##     ## ##     ## ##     ## ##     ## ########        
% ##     ## ##     ## ##     ## ##     ## ##   ##         
% ##     ## ##     ## ##     ## ##     ## ##    ##        
%  #######  ########   #######   #######  ##     ##       


% ########  ##     ## ##        ######  ########  ######  
% ##     ## ##     ## ##       ##    ## ##       ##    ## 
% ##     ## ##     ## ##       ##       ##       ##       
% ########  ##     ## ##        ######  ######    ######  
% ##        ##     ## ##             ## ##             ## 
% ##        ##     ## ##       ##    ## ##       ##    ## 
% ##         #######  ########  ######  ########  ######  

%%
% Can the DA Model account for responses to pulses of odor?





% ########     ##          ##     ##    ##     ##    ########     ######  
% ##     ##    ##          ##     ##    ###   ###    ##          ##    ## 
% ##     ##    ##          ##     ##    #### ####    ##          ##       
% ########     ##          ##     ##    ## ### ##    ######       ######  
% ##           ##          ##     ##    ##     ##    ##                ## 
% ##           ##          ##     ##    ##     ##    ##          ##    ## 
% ##           ########     #######     ##     ##    ########     ######                         