% fig_gain_control_fast_broadly_observed
% makes a figure for the paper showing that gain control is fast and is broadly observed


pHeader;

%% global parameters
history_lengths = round(logspace(1.7,4,50)); % all the history lengths we look at, in ms
example_history_length = 300; % this history length shown in the first row, in ms

figure('outerposition',[0 0 801 800],'PaperUnits','points','PaperSize',[801 800]); hold on


%  ######     ###    ########  ##        #######  ######## ########    ###    
% ##    ##   ## ##   ##     ## ##       ##     ##    ##       ##      ## ##   
% ##        ##   ##  ##     ## ##       ##     ##    ##       ##     ##   ##  
% ##       ##     ## ########  ##       ##     ##    ##       ##    ##     ## 
% ##       ######### ##   ##   ##       ##     ##    ##       ##    ######### 
% ##    ## ##     ## ##    ##  ##       ##     ##    ##       ##    ##     ## 
%  ######  ##     ## ##     ## ########  #######     ##       ##    ##     ## 

% ########     ###    ########    ###    
% ##     ##   ## ##      ##      ## ##   
% ##     ##  ##   ##     ##     ##   ##  
% ##     ## ##     ##    ##    ##     ## 
% ##     ## #########    ##    ######### 
% ##     ## ##     ##    ##    ##     ## 
% ########  ##     ##    ##    ##     ## 

clearvars -except being_published tau_h example_history_length

% load the data
if ~exist('orn_data','var')
	load('/local-data/DA-paper/fig4/Carlotta_Data.mat')
end



% ########  #### ######## ########         #######  ########   #######  ########   ######  
% ##     ##  ##  ##       ##              ##     ## ##     ## ##     ## ##     ## ##    ## 
% ##     ##  ##  ##       ##              ##     ## ##     ## ##     ## ##     ## ##       
% ##     ##  ##  ######   ######          ##     ## ##     ## ##     ## ########   ######  
% ##     ##  ##  ##       ##              ##     ## ##     ## ##     ## ##   ##         ## 
% ##     ##  ##  ##       ##       ###    ##     ## ##     ## ##     ## ##    ##  ##    ## 
% ########  #### ##       ##       ###     #######  ########   #######  ##     ##  ######  
do_these = [7 8 10 14 17];

ax(1) = subplot(2,2,1); hold on
ax(2) = subplot(2,2,3); hold on

c = lines(length(do_these));
for i = 1:length(do_these)
	temp = orn_data(do_these(i));
	pred = nanmean(temp.firing_projected,2); pred = pred(temp.use_this_segment);
	resp = nanmean(temp.firing_rate,2);  resp = resp(temp.use_this_segment);
	stim = nanmean(temp.stimulus,2); stim = stim(temp.use_this_segment);
 	stim = stim/nanmean(stim);

	% find when the valve opens
	[ons,offs] = findValveWhiffs(temp);

	% plot the gain in each of these windows
	[gain,gain_err] = findGainInWindows(ons,offs,pred,resp);

	gain = gain/nanmean(gain);

	rm_this = gain<0.4 | gain_err < .8;
	gain(rm_this) = [];
	ons(rm_this) = [];
	offs(rm_this) = [];

	% find the mean stimulus in the preceding X ms in these windows
	mean_stim = findMeanInWindows(ons,offs,computeSmoothedStimulus(stim,example_history_length));

	% % plot the gain vs. the mean stim after sorting it
	[mean_stim,idx] = sort(mean_stim);
	plot(ax(1),mean_stim,gain(idx),'+-','Color',c(i,:))

	% also find rho for various values of the history length and plot it
	rho = findRhoForHistoryLengths(gain,stim,ons,offs,history_lengths);

	plot(ax(2),history_lengths,rho,'.-','Color',c(i,:))
end
set(ax(2),'XScale','log','YLim',[-1 0])
set(ax(1),'XScale','log','YScale','log','XLim',[.1 2],'YLim',[.1 10])



lh = legend(h2,{'methyl butyrate','1-octen-3-ol','diethyl succinate','ethyl acetate','ethyl butyrate'});
lh.Position = [0.18    0.5081    0.0886    0.06];

% ########  #### ######## ########         #######  ########  ##    ##  ######  
% ##     ##  ##  ##       ##              ##     ## ##     ## ###   ## ##    ## 
% ##     ##  ##  ##       ##              ##     ## ##     ## ####  ## ##       
% ##     ##  ##  ######   ######          ##     ## ########  ## ## ##  ######  
% ##     ##  ##  ##       ##              ##     ## ##   ##   ##  ####       ## 
% ##     ##  ##  ##       ##       ###    ##     ## ##    ##  ##   ### ##    ## 
% ########  #### ##       ##       ###     #######  ##     ## ##    ##  ######  


do_these = [12 17];
example_orn = do_these(1);

% plot the responses to each valve on for one of them
[ons,offs] = findValveWhiffs(orn_data(example_orn));    % find times when valve opens
S = nanmean(orn_data(example_orn).stimulus,2);
use_this_segment = false(length(S),1);
for i = 1:length(ons)
	use_this_segment(ons(i):offs(i)) = true;
end
C = colourByMeanStimulus(S,use_this_segment,'history_length',example_history_length);  % build colour index based on mean stimulus
X = nanmean(orn_data(example_orn).firing_projected,2);
Y = nanmean(orn_data(example_orn).firing_rate,2);
X(~use_this_segment) = NaN;
Y(~use_this_segment) = NaN;
subplot(5,4,14), hold on 
scatter(X,Y,12,C,'filled')
xlabel('Projected Stimulus (V)')
ylabel('Firing Rate (Hz)')

ax(3) = subplot(5,4,15); hold on
ax(6) = subplot(5,4,16); hold on 

for i = 1:length(do_these)
	temp = orn_data(do_these(i));
	% now fit a NL
	% temp = fitNL(temp);
	plot(temp,[ax(3) ax(3) ax(6)],'valveGainAnalysis.firing_rate.mu','history_lengths',tau_h,'showNL',false,'history_length',example_history_length);
end

% also add the data with the naturalistic ab2 responses
temp = orn_data(25);
% temp = fitNL(temp);
plot(temp,[ax(3) ax(3) ax(6)],'excGainAnalysis.firing_rate.mu','history_lengths',round(logspace(2,4,30)),'showNL',false,'history_length',example_history_length,'normalise_gain',true);
% colour them nicely

h1 = get(ax(3),'Children');
h2 = get(ax(6),'Children');
% remove some extra random shit
delete(h1(2:end-2))
h1 = get(ax(3),'Children');
delete(h2(1));
h2 = get(ax(6),'Children');

c = lines(length(h1));

for i = 1:length(h2)
	set(h1(i),'Color',c(i,:),'MarkerFaceColor',c(i,:),'Marker','o')
	set(h2(i),'Color',c(i,:),'MarkerFaceColor',c(i,:),'Marker','o','LineStyle','-')
end
set(ax(3),'YScale','log','YLim',[.1 5])
xlabel(ax(3),'\mu_{stimulus} (norm)')
ylabel(ax(3),'Gain (norm)')
set(ax(6),'XScale','log','XLim',[10 1e4],'XTick',[1e1 1e2 1e3 1e4],'YLim',[-1 0])

lh = legend(h2,{'ab2A','ab3A','pb1A'});
lh.Position = [0.2    0.3517    0.0849    0.04];
uistack(h2(2),'top')

prettyFig('fs',12)
labelFigure

if being_published	
	snapnow	
	delete(gcf)
end

%% Supplementary Figure: 
% in this supplementary figure, we show the response vs. the projected stimulus, and colour the data by the mean stimulus in the preceding X ms

% plot the responses to each valve on for one of them
[ons,offs] = findValveWhiffs(orn_data(example_orn));    % find times when valve opens
S = nanmean(orn_data(example_orn).stimulus,2);
use_this_segment = false(length(S),1);
for i = 1:length(ons)
	use_this_segment(ons(i):offs(i)) = true;
end
C = colourByMeanStimulus(S,use_this_segment,'history_length',example_history_length);  % build colour index based on mean stimulus
X = nanmean(orn_data(example_orn).firing_projected,2);
Y = nanmean(orn_data(example_orn).firing_rate,2);
X(~use_this_segment) = NaN;
Y(~use_this_segment) = NaN;
subplot(5,4,14), hold on 
scatter(X,Y,12,C,'filled')
xlabel('Projected Stimulus (V)')
ylabel('Firing Rate (Hz)')

%% Supp Figure

% ########  ######## ########  ##       ####  ######     ###    ######## ########  ######  
% ##     ## ##       ##     ## ##        ##  ##    ##   ## ##      ##    ##       ##    ## 
% ##     ## ##       ##     ## ##        ##  ##        ##   ##     ##    ##       ##       
% ########  ######   ########  ##        ##  ##       ##     ##    ##    ######    ######  
% ##   ##   ##       ##        ##        ##  ##       #########    ##    ##             ## 
% ##    ##  ##       ##        ##        ##  ##    ## ##     ##    ##    ##       ##    ## 
% ##     ## ######## ##        ######## ####  ######  ##     ##    ##    ########  ######  


do_these = [7 9 13 15 16];
example_orn = do_these(2);

% plot the responses to each valve on for one of them
[ons,offs] = findValveWhiffs(orn_data(example_orn));    % find times when valve opens
S = nanmean(orn_data(example_orn).stimulus,2);
use_this_segment = false(length(S),1);
for i = 1:length(ons)
	use_this_segment(ons(i):offs(i)) = true;
end


for i = 1:length(do_these)
	temp = orn_data(do_these(i));
	plot(temp,[ax(1) ax(1) ax(4)],'valveGainAnalysis.firing_rate.mu','history_lengths',tau_h,'showNL',false,'history_length',example_history_length);
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
set(ax(4),'XScale','log','XLim',[10 1e4],'XTick',[1e1 1e2 1e3 1e4],'YLim',[-1 0])
xlabel(ax(1),'\mu_{stimulus} (norm)')
ylabel(ax(1),'Gain (norm)')

lh = legend(h2,{'5/28','6/05','6/12','6/19','6/19'});
lh.Position = [0.2   0.6807    0.0749    0.06];


%% Version Info
%
pFooter;



