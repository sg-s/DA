% DA_Paper.m
% 
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
%
% Outline: 
% 1) Do ORNs modulate gain?

% 1.1) Do ORNs modulate gain in flickering stimulus (1-octen-3-ol)
% 		We fit a linear filter --done
% 		Perform gain analysis vs. linear filter
% 		Bootstrap statistical significance 
% 		We conclude there is gain adaptation in this dataset

% 1.2) Do ORNs modulate gain in other flickering stim data? 
% 		fit liner filter (+recitfier) to all other data sets
% 		repeat gain analysis + Bootstrap statistics 



% 2) Can we explain fast gain adaptation? 

% 2.1) can a NLN model explain fast gain adaptation? 
% 		fit NLN model to all data (flickering, turbulent plume)
% 		build a table of parameters for all data
% 		do gain analysis + bootstrap statistics

% 2.2) can a DA model explain this fast gain adaptation? 
% 		fit DA model to all data (flickering, turbulent plume)
% 		build a table of parameters for all data
% 		gain analysis vs DA model + bootstrap statistics 
% 
% 3) How relevant is gain adaptation? Do ORNs modulate gain in turbulent air stream experiment? 
% 		fit linear model to all data 
% 		do gain analysis + bootstrap statistics for all data

%
%
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



% load data
load('/local-data/DA-paper/data.mat')

%% Do ORNs rapidly modulate gain?
% Pseudo-white-noise analysis of ORN responses involves presenting binary flickering pulses of odor to the ORN and recording their response. If ORNs rapidly modulate gain on the timescale of response, then responses to pulses of odor in the sequence where the stimulus is locally low will be different from responses to pulses of odor in the sequence where the stimulus is locally high. 

%%
% The data looks like this. The following figure shows the valve state, the odor concentration, and the neuron response. The neuron is ab3A, and the odor presented is 1-octen-3-ol diluted to 3x $10^{-3}$ in Paraffin Oil. The correlation time in the valve position is 30ms. 

%%
% This data file is used for the following analysis:
td = 7;

disp(data(td).original_name)

% detrend PID with a quadratic term
ptrend = fit(data(td).time(:),data(td).PID(:),'Poly2'); 
data(td).PID = data(td).PID(:) - (ptrend(data(td).time) - mean(ptrend(data(td).time)));


% build a simple linear model
[K,~,filtertime] = FindBestFilter(data(td).PID(500:end),data(td).ORN(500:end),[],'filter_length=201;');
data(td).K = K;
data(td).filtertime = filtertime*mean(diff(data(td).time));
data(td).LinearFit = mean(data(td).ORN) + convolve(data(td).time,data(td).PID,data(td).K,data(td).filtertime);




figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(2,1,1), hold on
plot(data(td).time,data(td).PID,'k');
set(gca,'XLim',[min(data(td).time) max(data(td).time)])
ylabel('PID Voltage (V)')

subplot(2,1,2), hold on
plot(data(td).time,data(td).ORN,'k');
plot(data(td).time,data(td).LinearFit,'r');
set(gca,'XLim',[min(data(td).time) max(data(td).time)])
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
legend ORN LinearFit
PrettyFig;

figure('outerposition',[0 0 400 400],'PaperUnits','points','PaperSize',[1000 400]); hold on
plot(data(td).filtertime,data(td).K,'k','LineWidth',2)
title('Linear Filter for this data')
xlabel('FitlerLag (s)')
set(gca,'XLim',[min(data(td).filtertime) max(data(td).filtertime)])
PrettyFig;


%%
% The distributions of the input to neuron and the neuron responses are shown below on the left. The autocorrelaiton functions of the input and output are shown on the right. If the ORN modulates its gain on a rapid time-scale, it must do so in this case on a time-scale smaller than the autocorrelation time of the stimulus. 

[act] = PlotDataStatistics(data,td);


%%
% Does the ORN selectively amplify responses to relatively small stimuli and suppress responses to relatively high stimuli?
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ph(3) = subplot(1,2,1); hold on 
axis square
ph(4) = subplot(1,2,2); hold on

s = 300; % when we start for the gain analysis
z = length(data(td).ORN); % where we end
example_history_length = [0.09];
history_lengths = [.030:0.06:2*act];

clear x
x.response = data(td).ORN(s:z);
x.prediction = data(td).LinearFit(s:z);
x.stimulus = data(td).PID(s:z);
x.time = data(td).time(s:z);
x.filter_length = 201;

redo_bootstrap = 0;
if redo_bootstrap
	data(td).LinearFit_p = GainAnalysis3(x,history_lengths,example_history_length,ph);
else
	GainAnalysis3(x,history_lengths,example_history_length,ph,data(td).LinearFit_p);
end
clear x

% save
save('/local-data/DA-paper/data.mat','data','-append')

%%
% We now repeat the analysis on a different data set, where the same odor is presented to the same type of neuron, but the flickering stimulus has a longer correlation time. The response of the neuron is now quite different, and sometimes goes to 0 (stops firing entirely). 

td = 9;
%%
% This data file is being used for this analysis:
disp(data(td).original_name)



% detrend PID with a quadratic term
ptrend = fit(data(td).time(:),data(td).PID(:),'Poly2'); 
data(td).PID = data(td).PID(:) - (ptrend(data(td).time) - mean(ptrend(data(td).time)));


% build a simple linear model
[K,~,filtertime] = FindBestFilter(data(td).PID(500:end),data(td).ORN(500:end),[],'filter_length=201;','min_cutoff = 0;');
data(td).K = K;
data(td).filtertime = filtertime*mean(diff(data(td).time));
data(td).LinearFit = mean(data(td).ORN) + convolve(data(td).time,data(td).PID,data(td).K,data(td).filtertime);
data(td).LinearFit(data(td).LinearFit <0 ) = 0;


%%
% This is what the data looks like. The top panel shows the stimulus and the bottom panel shows the neuron response together with the linear fit. 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(2,1,1), hold on
plot(data(td).time,data(td).PID,'k');
set(gca,'XLim',[min(data(td).time) max(data(td).time)])
ylabel('PID Voltage (V)')

subplot(2,1,2), hold on
plot(data(td).time,data(td).ORN,'k');
plot(data(td).time,data(td).LinearFit,'r');
set(gca,'XLim',[min(data(td).time) max(data(td).time)])
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
legend ORN LinearFit
PrettyFig;

%%
% How is the filter changed? The following figure shows the filter computed from this dataset compared to previous dataset.
figure('outerposition',[0 0 400 400],'PaperUnits','points','PaperSize',[1000 400]); hold on
plot(data(td).filtertime,data(td).K,'r','LineWidth',2), hold on
plot(data(td).filtertime,data(7).K,'k','LineWidth',2)
title('Linear Filter for this data')
xlabel('FitlerLag (s)')
set(gca,'XLim',[min(data(td).filtertime) max(data(td).filtertime)])
PrettyFig;
legend ThisData PreviousData


% plot statistics of the data: histograms and auto-correlation functions
[act] = PlotDataStatistics(data,td);


% We can perform a similar gain analysis like we did on the synthetic data on the real data from the ORN. The plot on the left compares the ORN response data (on the Y-axis) to the linear fit, while the plot on the right compares the data to the DA model fit. 
s = 300; % when we start for the gain analysis
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ph(3) = subplot(1,2,1); hold on 
axis square
ph(4) = subplot(1,2,2); hold on
z = length(data(td).ORN); % where we end
example_history_length = [0.3900];
history_lengths = [.030:0.06:2*act];

clear x
x.response = data(td).ORN(s:z);
x.prediction = data(td).LinearFit(s:z);
x.stimulus = data(td).PID(s:z);
x.time = data(td).time(s:z);
x.filter_length = 201;

redo_bootstrap = 0;
if redo_bootstrap
	data(td).LinearFit_p = GainAnalysis3(x,history_lengths,example_history_length,ph);
else
	if isempty(data(td).LinearFit_p)
		data(td).LinearFit_p = GainAnalysis3(x,history_lengths,example_history_length,ph);
	else
		GainAnalysis3(x,history_lengths,example_history_length,ph,data(td).LinearFit_p);
	end
end
clear x

% save
save('/local-data/DA-paper/data.mat','data','-append')

return






%Now that we have the filter, fit it to dose-response data
data_new=load('/local-data/orn/carlotta-paper/fig3abc.mat');
data(2).stimulus = data_new.stimulus;
data(2).time = data_new.time;
data(2).response = data_new.response;
data(2).DAFit = data_new.DAFit;
data(2).p = data_new.p;
clear data_new

% convert filter into new timescale
nt = min(filtertime):mean(diff(data(2).time)):max(filtertime);
data(2).K = interp1(filtertime,K,nt);

% fit linear model
[data(2).LinearFit data(2).K] = FitLinearModel2Data(data(2),data(2).K);
data(2).LinearFit = reshape(data(2).LinearFit,length(data(2).response),10);

% make the plot
figure, hold on
subplot(1,4,1), hold on
plot(data(2).time,data(2).stimulus)
subplot(1,4,2), hold on
plot(data(2).time,data(2).response)
subplot(1,4,3), hold on
plot(data(2).time,data(2).LinearFit)
subplot(1,4,4), hold on
plot(data(2).time,data(2).DAFit)



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

%% Can the DA Model account for fast adaptation? 
% Here we fit the DA model to the data. 


return

if ~isfield(data(1),'DAFit')
	x0 = [5214  3  0.35 1.8 2   9   2  20];
	lb = [1200  0    0    0   2   0   2  0];
	ub = [9200  300  1    10  2  30   2  30];
	psoptions = psoptimset('UseParallel',true,'CompletePoll', 'on', 'Vectorized', 'off','Display','iter','MaxIter',2000,'MaxFunEvals',10000);
	x = patternsearch(@(x) DA_cost_function(x,data(1),@Cost,500,1),x0,[],[],[],[],lb,ub,psoptions);
	data(1).DAFit_p = ValidateDAParameters2(x);
	pred = DA_integrate2(data(1).stimulus,data(1).DAFit_p);
	pred = pred-mean(pred);
	multiplot([],data(1).response,pred)
	clear x
	data(1).DAFit_p
	% save('DA_Paper_data.mat','data','-append')
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