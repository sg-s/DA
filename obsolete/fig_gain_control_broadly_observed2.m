% fast gain control is broadly observed 
% makes figure 6 of the paper, showing that gain control is broadly observed
% 
% created by Srinivas Gorur-Shandilya at 9:56 , 04 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

pHeader;

% load the data
if ~exist('orn_data','var')
	load('Carlotta_Data.mat')
end

% % detrend the stimulus for carlotta's data
% for i = 1:length(orn_data)
% 	if ~isempty(strfind(orn_data(i).data_creator,'Carlotta'))
% 		S = orn_data(i).stimulus;
% 		for j = 1:width(S)
% 			ff = fit((1:length(S))',S(:,j),'poly3');
% 			S(:,j) = S(:,j) - ff(1:length(S)) + mean(S(:,j));
% 		end
% 		orn_data(i).stimulus = S;
% 	end
% end

history_lengths = round(logspace(1.7,4,30));

figure('outerposition',[0 0 1400 800],'PaperUnits','points','PaperSize',[1400 800]); hold on
clear ax
for i = 1:6
	ax(i) = subplot(2,3,i); hold on
end

% ########  ######## ########  ##       ####  ######     ###    ######## ########  ######  
% ##     ## ##       ##     ## ##        ##  ##    ##   ## ##      ##    ##       ##    ## 
% ##     ## ##       ##     ## ##        ##  ##        ##   ##     ##    ##       ##       
% ########  ######   ########  ##        ##  ##       ##     ##    ##    ######    ######  
% ##   ##   ##       ##        ##        ##  ##       #########    ##    ##             ## 
% ##    ##  ##       ##        ##        ##  ##    ## ##     ##    ##    ##       ##    ## 
% ##     ## ######## ##        ######## ####  ######  ##     ##    ##    ########  ######  


do_these = [7 9 13 15 16];
for i = 1:length(do_these)
	% now fit a NL
	temp = orn_data(do_these(i));
	temp = fitNL(temp);
	plot(temp,[ax(1) ax(1) ax(4)],'valveGainAnalysis.firing_rate.mu','history_lengths',history_lengths,'showNL',false,'history_length',300);
end
% colour them nicely
c = lines(length(do_these));
h1 = get(ax(1),'Children');
h2 = get(ax(4),'Children');
for i = 1:length(do_these)
	set(h1(i),'Color',c(i,:),'MarkerFaceColor',c(i,:),'Marker','o')
	set(h2(i),'Color',c(i,:),'MarkerFaceColor',c(i,:),'Marker','o','LineStyle','-')
end
set(ax(1),'YScale','log','YLim',[.1 5])
xlabel(ax(1),'\mu_{stimulus} (norm)')
ylabel(ax(1),'Gain (norm)')
title(ax(1),'Experimental Replicates')

legend(h2,{'5/28','6/05','6/12','6/19','6/19'})

% ########  #### ######## ########         #######  ########   #######  ########   ######  
% ##     ##  ##  ##       ##              ##     ## ##     ## ##     ## ##     ## ##    ## 
% ##     ##  ##  ##       ##              ##     ## ##     ## ##     ## ##     ## ##       
% ##     ##  ##  ######   ######          ##     ## ##     ## ##     ## ########   ######  
% ##     ##  ##  ##       ##              ##     ## ##     ## ##     ## ##   ##         ## 
% ##     ##  ##  ##       ##       ###    ##     ## ##     ## ##     ## ##    ##  ##    ## 
% ########  #### ##       ##       ###     #######  ########   #######  ##     ##  ######  

do_these = [7 8 10 14 17];
for i = 1:length(do_these)
	% now fit a NL
	temp = orn_data(do_these(i));
	temp = fitNL(temp);
	plot(temp,[ax(2) ax(2) ax(5)],'valveGainAnalysis.firing_rate.mu','history_lengths',history_lengths,'showNL',false,'history_length',300);
end
% colour them nicely
c = lines(length(do_these));
h1 = get(ax(2),'Children');
h2 = get(ax(5),'Children');
for i = 1:length(h2)
	set(h1(i),'Color',c(i,:),'MarkerFaceColor',c(i,:),'Marker','o')
	set(h2(i),'Color',c(i,:),'MarkerFaceColor',c(i,:),'Marker','o','LineStyle','-')
end
set(ax(2),'YScale','log','YLim',[.1 5])
xlabel(ax(2),'\mu_{stimulus} (norm)')
ylabel(ax(2),'Gain (norm)')
title(ax(2),'Diff. odors')

legend(h2,{'1but','1o3ol','d2succ','2ac','2but'})


% ########  #### ######## ########         #######  ########  ##    ##  ######  
% ##     ##  ##  ##       ##              ##     ## ##     ## ###   ## ##    ## 
% ##     ##  ##  ##       ##              ##     ## ##     ## ####  ## ##       
% ##     ##  ##  ######   ######          ##     ## ########  ## ## ##  ######  
% ##     ##  ##  ##       ##              ##     ## ##   ##   ##  ####       ## 
% ##     ##  ##  ##       ##       ###    ##     ## ##    ##  ##   ### ##    ## 
% ########  #### ##       ##       ###     #######  ##     ## ##    ##  ######  

do_these = [12 17];
for i = 1:length(do_these)
	% now fit a NL
	temp = orn_data(do_these(i));
	temp = fitNL(temp);
	plot(temp,[ax(3) ax(3) ax(6)],'valveGainAnalysis.firing_rate.mu','history_lengths',history_lengths,'showNL',false,'history_length',300);
end
temp = orn_data(25);
temp = fitNL(temp);
plot(temp,[ax(3) ax(3) ax(6)],'excGainAnalysis.firing_rate.mu','history_lengths',round(logspace(2,4,30)),'showNL',false,'history_length',300);
% colour them nicely

h1 = get(ax(3),'Children');
c = lines(length(h1));
h2 = get(ax(6),'Children');
for i = 1:length(h2)
	set(h1(i),'Color',c(i,:),'MarkerFaceColor',c(i,:),'Marker','o')
	set(h2(i),'Color',c(i,:),'MarkerFaceColor',c(i,:),'Marker','o','LineStyle','-')
end
set(ax(3),'YScale','log','YLim',[.1 5])
xlabel(ax(3),'\mu_{stimulus} (norm)')
ylabel(ax(3),'Gain (norm)')
title(ax(3),'Diff. ORNs')

legend(h2,{'ab2A','ab3A','pb1A'})

prettyFig('plw=1.5;','lw=1.5;','fs=14;')


if being_published
	snapnow
	delete(gcf)
end
return

%% -ve Controls: LN Model
% In this section, we generate synthetic data using the LN model, and re-run this analysis to see if we pick something up. 

if ~exist('LN_data','var')
	LN_data = ORNData;
	for i = 1:length(orn_data)
		textbar(i,length(orn_data))
		temp = orn_data(i);
		temp = fitNL(temp);

		LN_data(i).stimulus = nanmean(temp.stimulus,2);
		LN_data(i).firing_rate = nanmean(temp.firing_projected,2);
		LN_data(i).valve = temp.valve;
		LN_data(i).K_firing = nanmean(temp.K_firing,2);
		LN_data(i).filtertime_firing = temp.filtertime_firing;
	end

	% back out filters everywhere
	for i = 1:length(LN_data)
		disp(i)
		LN_data(i).regularisation_factor = 1;
		LN_data(i) = backOutFilters(LN_data(i));
	end

end


figure('outerposition',[0 0 1400 800],'PaperUnits','points','PaperSize',[1400 800]); hold on
clear ax
for i = 1:6
	ax(i) = subplot(2,3,i); hold on
end

% ########  ######## ########  ##       ####  ######     ###    ######## ########  ######  
% ##     ## ##       ##     ## ##        ##  ##    ##   ## ##      ##    ##       ##    ## 
% ##     ## ##       ##     ## ##        ##  ##        ##   ##     ##    ##       ##       
% ########  ######   ########  ##        ##  ##       ##     ##    ##    ######    ######  
% ##   ##   ##       ##        ##        ##  ##       #########    ##    ##             ## 
% ##    ##  ##       ##        ##        ##  ##    ## ##     ##    ##    ##       ##    ## 
% ##     ## ######## ##        ######## ####  ######  ##     ##    ##    ########  ######  


do_these = [7 9 13 15 16];
for i = 1:length(do_these)
	% now fit a NL
	temp = LN_data(do_these(i));
	temp = fitNL(temp);
	plot(temp,[ax(1) ax(1) ax(4)],'valveGainAnalysis.firing_rate.mu','history_lengths',history_lengths,'showNL',false,'history_length',300);
end
% colour them nicely
c = lines(length(do_these));
h1 = get(ax(1),'Children');
h2 = get(ax(4),'Children');
for i = 1:length(do_these)
	set(h1(i),'Color',c(i,:),'MarkerFaceColor',c(i,:),'Marker','o')
	set(h2(i),'Color',c(i,:),'MarkerFaceColor',c(i,:),'Marker','o','LineStyle','-')
end
set(ax(1),'YScale','log','YLim',[.1 5])
xlabel(ax(1),'\mu_{stimulus} (norm)')
ylabel(ax(1),'Gain (norm)')
title(ax(1),'Experimental Replicates')

legend(h2,{'5/28','6/05','6/12','6/19','6/19'})


% ########  #### ######## ########         #######  ########   #######  ########   ######  
% ##     ##  ##  ##       ##              ##     ## ##     ## ##     ## ##     ## ##    ## 
% ##     ##  ##  ##       ##              ##     ## ##     ## ##     ## ##     ## ##       
% ##     ##  ##  ######   ######          ##     ## ##     ## ##     ## ########   ######  
% ##     ##  ##  ##       ##              ##     ## ##     ## ##     ## ##   ##         ## 
% ##     ##  ##  ##       ##       ###    ##     ## ##     ## ##     ## ##    ##  ##    ## 
% ########  #### ##       ##       ###     #######  ########   #######  ##     ##  ######  

do_these = [7 8 10 14 17];
for i = 1:length(do_these)
	% now fit a NL
	temp = LN_data(do_these(i));
	temp = fitNL(temp);
	plot(temp,[ax(2) ax(2) ax(5)],'valveGainAnalysis.firing_rate.mu','history_lengths',history_lengths,'showNL',false,'history_length',300);
end
% colour them nicely
c = lines(length(do_these));
h1 = get(ax(2),'Children');
h2 = get(ax(5),'Children');
for i = 1:length(h2)
	set(h1(i),'Color',c(i,:),'MarkerFaceColor',c(i,:),'Marker','o')
	set(h2(i),'Color',c(i,:),'MarkerFaceColor',c(i,:),'Marker','o','LineStyle','-')
end
set(ax(2),'YScale','log','YLim',[.1 5])
xlabel(ax(2),'\mu_{stimulus} (norm)')
ylabel(ax(2),'Gain (norm)')
title(ax(2),'Diff. odors')

legend(h2,{'1but','1o3ol','d2succ','2ac','2but'})



% ########  #### ######## ########         #######  ########  ##    ##  ######  
% ##     ##  ##  ##       ##              ##     ## ##     ## ###   ## ##    ## 
% ##     ##  ##  ##       ##              ##     ## ##     ## ####  ## ##       
% ##     ##  ##  ######   ######          ##     ## ########  ## ## ##  ######  
% ##     ##  ##  ##       ##              ##     ## ##   ##   ##  ####       ## 
% ##     ##  ##  ##       ##       ###    ##     ## ##    ##  ##   ### ##    ## 
% ########  #### ##       ##       ###     #######  ##     ## ##    ##  ######  

do_these = [12 17];
for i = 1:length(do_these)
	% now fit a NL
	temp = LN_data(do_these(i));
	temp = fitNL(temp);
	plot(temp,[ax(3) ax(3) ax(6)],'valveGainAnalysis.firing_rate.mu','history_lengths',history_lengths,'showNL',false,'history_length',300);
end
temp = LN_data(25);
temp = fitNL(temp);
plot(temp,[ax(3) ax(3) ax(6)],'excGainAnalysis.firing_rate.mu','showNL',false,'history_length',300);
% colour them nicely

h1 = get(ax(3),'Children');
c = lines(length(h1));
h2 = get(ax(6),'Children');
for i = 1:length(h2)
	set(h1(i),'Color',c(i,:),'MarkerFaceColor',c(i,:),'Marker','o')
	set(h2(i),'Color',c(i,:),'MarkerFaceColor',c(i,:),'Marker','o','LineStyle','-')
end
set(ax(3),'YScale','log','YLim',[.1 5])
xlabel(ax(3),'\mu_{stimulus} (norm)')
ylabel(ax(3),'Gain (norm)')
title(ax(3),'Diff. ORNs')

legend(h2,{'ab2A','ab3A','pb1A'})

prettyFig('plw=1.5;','lw=1.5;','fs=14;')


if being_published
	snapnow
	delete(gcf)
end


%% +ve Controls: DA Model
% In this section, we generate synthetic data using the DA model, and re-run this analysis to see if we pick something up. 

clear p0
p0.   s0 = -0.0617;
p0.  n_z = 2;
p0.tau_z = 55.1249;
p0.  n_y = 2.7500;
p0.tau_y = 10.7002;
p0.    C = 0.0144;
p0.    A = 229.2877;
p0.    B = 3.6578;

if ~exist('DA_data','var')

	% fit DA models to each dataset
	% for i = 1:length(orn_data)
	% 	clear d
	% 	d.response = nanmean(orn_data(i).firing_rate,2);
	% 	d.response(1:1e3) = NaN;
	% 	S = nanmean(orn_data(i).stimulus,2);
	% 	S = S - min(S); S = S/max(S);
	% 	d.stimulus = S;
	% 	p(i) = fitModel2Data(@DAModelv2,d,'p0',p0);
	% end
	load('.cache/carlotta_DA_model_fits.mat','p')

	DA_data = ORNData;
	for i = 1:25
		textbar(i,length(orn_data))
		temp = orn_data(i);
		S = nanmean(temp.stimulus,2);
		S = S - min(S); S = S/max(S);
		DA_data(i).stimulus = S;
		R = DAModelv2(DA_data(i).stimulus,p(i)); R(R<0) = 0;
		DA_data(i).firing_rate = R;
		DA_data(i).valve = temp.valve;
	end

	% back out filters everywhere
	for i = 1:25
		disp(i)
		DA_data(i).regularisation_factor = 1;
		DA_data(i) = backOutFilters(DA_data(i));
	end

end


figure('outerposition',[0 0 1400 800],'PaperUnits','points','PaperSize',[1400 800]); hold on
clear ax
for i = 1:6
	ax(i) = subplot(2,3,i); hold on
end

% ########  ######## ########  ##       ####  ######     ###    ######## ########  ######  
% ##     ## ##       ##     ## ##        ##  ##    ##   ## ##      ##    ##       ##    ## 
% ##     ## ##       ##     ## ##        ##  ##        ##   ##     ##    ##       ##       
% ########  ######   ########  ##        ##  ##       ##     ##    ##    ######    ######  
% ##   ##   ##       ##        ##        ##  ##       #########    ##    ##             ## 
% ##    ##  ##       ##        ##        ##  ##    ## ##     ##    ##    ##       ##    ## 
% ##     ## ######## ##        ######## ####  ######  ##     ##    ##    ########  ######  


do_these = [7 9 13 15 16];
for i = 1:length(do_these)
	% now fit a NL
	temp = DA_data(do_these(i));
	temp = fitNL(temp);
	plot(temp,[ax(1) ax(1) ax(4)],'valveGainAnalysis.firing_rate.mu','history_lengths',history_lengths,'showNL',false,'history_length',300);
end
% colour them nicely
c = lines(length(do_these));
h1 = get(ax(1),'Children');
h2 = get(ax(4),'Children');
for i = 1:length(do_these)
	set(h1(i),'Color',c(i,:),'MarkerFaceColor',c(i,:),'Marker','o')
	set(h2(i),'Color',c(i,:),'MarkerFaceColor',c(i,:),'Marker','o','LineStyle','-')
end
set(ax(1),'YScale','log','YLim',[.1 5])
xlabel(ax(1),'\mu_{stimulus} (norm)')
ylabel(ax(1),'Gain (norm)')
title(ax(1),'Experimental Replicates')

legend(h2,{'5/28','6/05','6/12','6/19','6/19'})


% ########  #### ######## ########         #######  ########   #######  ########   ######  
% ##     ##  ##  ##       ##              ##     ## ##     ## ##     ## ##     ## ##    ## 
% ##     ##  ##  ##       ##              ##     ## ##     ## ##     ## ##     ## ##       
% ##     ##  ##  ######   ######          ##     ## ##     ## ##     ## ########   ######  
% ##     ##  ##  ##       ##              ##     ## ##     ## ##     ## ##   ##         ## 
% ##     ##  ##  ##       ##       ###    ##     ## ##     ## ##     ## ##    ##  ##    ## 
% ########  #### ##       ##       ###     #######  ########   #######  ##     ##  ######  

do_these = [7 8 10 14 17];
for i = 1:length(do_these)
	% now fit a NL
	temp = DA_data(do_these(i));
	temp = fitNL(temp);
	plot(temp,[ax(2) ax(2) ax(5)],'valveGainAnalysis.firing_rate.mu','history_lengths',history_lengths,'showNL',false,'history_length',300);
end
% colour them nicely
c = lines(length(do_these));
h1 = get(ax(2),'Children');
h2 = get(ax(5),'Children');
for i = 1:length(h2)
	set(h1(i),'Color',c(i,:),'MarkerFaceColor',c(i,:),'Marker','o')
	set(h2(i),'Color',c(i,:),'MarkerFaceColor',c(i,:),'Marker','o','LineStyle','-')
end
set(ax(2),'YScale','log','YLim',[.1 5])
xlabel(ax(2),'\mu_{stimulus} (norm)')
ylabel(ax(2),'Gain (norm)')
title(ax(2),'Diff. odors')

legend(h2,{'1but','1o3ol','d2succ','2ac','2but'})



% ########  #### ######## ########         #######  ########  ##    ##  ######  
% ##     ##  ##  ##       ##              ##     ## ##     ## ###   ## ##    ## 
% ##     ##  ##  ##       ##              ##     ## ##     ## ####  ## ##       
% ##     ##  ##  ######   ######          ##     ## ########  ## ## ##  ######  
% ##     ##  ##  ##       ##              ##     ## ##   ##   ##  ####       ## 
% ##     ##  ##  ##       ##       ###    ##     ## ##    ##  ##   ### ##    ## 
% ########  #### ##       ##       ###     #######  ##     ## ##    ##  ######  

do_these = [12 17];
for i = 1:length(do_these)
	% now fit a NL
	temp = DA_data(do_these(i));
	temp = fitNL(temp);
	plot(temp,[ax(3) ax(3) ax(6)],'valveGainAnalysis.firing_rate.mu','history_lengths',history_lengths,'showNL',false,'history_length',300);
end
temp = DA_data(25);
temp = fitNL(temp);
plot(temp,[ax(3) ax(3) ax(6)],'excGainAnalysis.firing_rate.mu','history_lengths',round(logspace(2,4,30)),'showNL',false,'history_length',300);
% colour them nicely

h1 = get(ax(3),'Children');
c = lines(length(h1));
h2 = get(ax(6),'Children');
for i = 1:length(h2)
	set(h1(i),'Color',c(i,:),'MarkerFaceColor',c(i,:),'Marker','o')
	set(h2(i),'Color',c(i,:),'MarkerFaceColor',c(i,:),'Marker','o','LineStyle','-')
end
set(ax(3),'YScale','log','YLim',[.1 5])
xlabel(ax(3),'\mu_{stimulus} (norm)')
ylabel(ax(3),'Gain (norm)')
title(ax(3),'Diff. ORNs')

legend(h2,{'ab2A','ab3A','pb1A'})

prettyFig('plw=1.5;','lw=1.5;','fs=14;')


if being_published
	snapnow
	delete(gcf)
end

% ##        #######   ######  ##     ##  ######  ######## 
% ##       ##     ## ##    ## ##     ## ##    ##    ##    
% ##       ##     ## ##       ##     ## ##          ##    
% ##       ##     ## ##       ##     ##  ######     ##    
% ##       ##     ## ##       ##     ##       ##    ##    
% ##       ##     ## ##    ## ##     ## ##    ##    ##    
% ########  #######   ######   #######   ######     ##    

% % load data
% if ~exist('locust_data','var')
% 	load('/local-data/DA-paper/locust/example-data')

% 	% clean up, sub-sample to 1ms
% 	PID = PID1; clear PID1
% 	EAG = EAG1; clear EAG1 

% 	PID = PID(:,1:10:end)';
% 	EAG = EAG(:,1:10:end)';
% 	valve = ODR1(:,1:10:end)';
% 	valve(valve<max(max(valve))/2) = 0;
% 	valve(valve>0) = 1;

% 	% set zero
% 	for i = 1:width(PID)
% 		PID(:,i) = PID(:,i) - mean(PID(1:300,i));
% 		EAG(:,i) = EAG(:,i) - mean(EAG(1:300,i));
% 		% filter
% 		PID(:,i) = bandPass(PID(:,i),Inf,30);
% 		EAG(:,i) = bandPass(EAG(:,i),2e3,Inf);
% 	end


% 	locust_data = ORNData;
% 	locust_data.filtertime_LFP = -.5:1e-3:.9;
% 	locust_data.regularisation_factor = 1;
% 	locust_data.stimulus = PID;
% 	locust_data.LFP = EAG;
% 	locust_data = backOutFilters(locust_data);
% end



%% Version Info
% The file that generated this document is called:
pFooter;
