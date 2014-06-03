% DA_Paper.m
% How this is organised:
% all the data (every single bit in every figure), is in a structure called "data"
% data(1) is response of ab3A to a flickering pulse sequence of 1-octen-3-ol 






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
% Pseudo-white-noise analysis of ORN responses involves presenting binary flickering pulses of odor to the ORN and recording their response. If ORNs rapidly modulate gain on the timescale of response, then responses to pulses of odor in the sequence where the stimulus is locally low will be different from responses to pulses of odor in the sequence where the stimulus is locally high. 

%%
% The data looks like this. The following figure shows the valve state, the odor concentration, and the neuron response. The neuron is ab3A, and the odor presented is 1-octen-3-ol diluted to 3x10^-^3 in Paraffin Oil. The correlation time in the valve position is 30ms. 


if ~exist('PID','var')
	filename = '/local-data/DA-paper/fig1/final_2011_06_14_ab3A_1o3ol3X-3_20ml_30sec_30ms_rand.mat';
	[PID, time, f,Valve,uncropped] = PrepData3(filename);
	PID = PID(:);
	time = time(:);
	f = f(:);

	% detrend PID
	ptrend = fit(time,PID,'Poly1'); 
	PID = PID - (ptrend(time) - mean(ptrend(time)));


	% assemble into a data structure
	data(1).stimulus = PID;
	data(1).response = f;
	data(1).Valve = Valve;
	data(1).time = time;
end


% build a simple linear model
[K,~,filtertime] = FindBestFilter(PID(500:end),f(500:end),[],'filter_length=201;');
data(1).K = K;
LinearFit = filter(K,1,PID-mean(PID))+ mean(f(500:end));
LinearFit(LinearFit<min(f)) = min(f);
% correct for offset
offset = abs(min(filtertime));
LinearFit(1:offset) = [];
LinearFit = [LinearFit; NaN(offset,1)];

% add it to the data
data(1).LinearFit = LinearFit;

figure('outerposition',[0 0 1500 1000],'PaperUnits','points','PaperSize',[1000 500]); hold on
fig(1).a(1) = subplot(3,6,2:6);
fig(1).a(2) = subplot(3,6,8:12);
fig(1).a(3) = subplot(3,6,13);
fig(1).a(4) = subplot(3,6,14:18); hold on

filtertime = filtertime*mean(diff(time));
offset = offset*mean(diff(time));

plot(fig(1).a(1),time,Valve);
plot(fig(1).a(2),time,PID);
plot(fig(1).a(3),filtertime,K);
plot(fig(1).a(4),time,f,'k'); 
plot(fig(1).a(4),time,LinearFit,'r')
set(fig(1).a(4),'XLim',[20 25])
set(fig(1).a(1),'XLim',[20 25],'YLim',[-0.1 1.1])
set(fig(1).a(2),'XLim',[20 25])
set(fig(1).a(3),'XLim',[min(filtertime-offset) max(filtertime)+offset])
PrettyFig;


% Real Data: Gain Analysis
% We can perform a similar gain analysis like we did on the synthetic data on the real data from the ORN. The plot on the left compares the ORN response data (on the Y-axis) to the linear fit, while the plot on the right compares the data to the DA model fit. 
s = 300; % when we start for the gain analysis
z = length(f); % where we end
example_history_length = [0.09];
history_lengths = [.030:0.06:1.002];

x.response = f(s:z);
x.prediction = LinearFit(s:z);
x.stimulus = PID(s:z);
x.time = time(s:z);
x.filter_length = 201;

GainAnalysis3(x,history_lengths,example_history_length);


return

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