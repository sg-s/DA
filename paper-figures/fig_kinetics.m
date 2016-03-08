% fig_kinetics.m
% 
% created by Srinivas Gorur-Shandilya at 2:14 , 07 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


pHeader;

figure('outerposition',[0 0 1100 700],'PaperUnits','points','PaperSize',[1100 700]); hold on
for i = 1:6
	ax(i) = subplot(2,3,i); hold on
end
delete(ax(4))
ax(1) = subplot(1,3,1); hold on


% show LFP slowdown 

[PID, LFP, fA, paradigm, orn, fly, AllControlParadigms, paradigm_hashes] = consolidateData('/local-data/DA-paper/LFP-MSG/september',1);

% remove baseline from all PIDs
for i = 1:width(PID)
	PID(:,i) = PID(:,i) - mean(PID(1:5e3,i));
end

% sort the paradigms sensibly
sort_value = [];
for i = 1:length(AllControlParadigms)
	sort_value(i) = (mean(AllControlParadigms(i).Outputs(1,:)));
end
[~,idx] = sort(sort_value);


AllControlParadigms = AllControlParadigms(idx);
paradigm_new = paradigm*NaN;
for i = 1:length(idx)
	paradigm_new(paradigm == idx(i)) = i;
end
paradigm = paradigm_new;

% throw our bad traces
bad_trials = (sum(fA) == 0 | isnan(sum(fA)));
PID(:,bad_trials) = [];
LFP(:,bad_trials) = [];
fA(:,bad_trials) = [];
paradigm(bad_trials) = [];
orn(bad_trials) = [];

% band pass all the LFP
filtered_LFP = LFP;
for i = 1:width(LFP)
	filtered_LFP(:,i) = filtered_LFP(:,i) - fastFiltFilt(ones(1e4,1),1e4,filtered_LFP(:,i));
	filtered_LFP(:,i) = filtered_LFP(:,i)*10; % to get the units right, now in mV
end

% some core variables
dt = 1e-3;
c = parula(max(paradigm)+1); % colour scheme


x_odor_LFP = NaN(1,width(PID));
x_odor_fA = NaN(1,width(PID));
x_LFP_fA = NaN(1,width(PID));

for i = 1:width(PID)
	s = PID(25e3:45e3,i)-mean(PID(25e3:45e3,i)); s = s/std(s);
	f = fA(25e3:45e3,i)-mean(fA(25e3:45e3,i)); f = f/std(f);
	r = filtered_LFP(25e3:45e3,i)-mean(filtered_LFP(25e3:45e3,i)); r = r/std(r); r = -r;

	temp = xcorr(f,s);
	[~,x_odor_fA(i)] = max(temp(19.5e3+1:20.5e3));

	temp = xcorr(r,s);
	[~,x_odor_LFP(i)] = max(temp(19.5e3+1:20.5e3));

	temp = xcorr(f,r);
	[~,x_LFP_fA(i)] = max(temp(19.5e3+1:20.5e3));
end

x_odor_LFP = x_odor_LFP - 500;
x_odor_fA = x_odor_fA - 500;
x_LFP_fA = x_LFP_fA - 500;

x_LFP_fA(x_LFP_fA<-100) = NaN;
x_odor_fA(x_odor_fA<-100) = NaN;
x_odor_LFP(x_odor_LFP<-100) = NaN;
x_LFP_fA(x_LFP_fA>0) = NaN;

%  ######  ######## #### ##     ##         ##    ##       ######## ########  
% ##    ##    ##     ##  ###   ###          ##   ##       ##       ##     ## 
% ##          ##     ##  #### ####           ##  ##       ##       ##     ## 
%  ######     ##     ##  ## ### ## #######    ## ##       ######   ########  
%       ##    ##     ##  ##     ##           ##  ##       ##       ##        
% ##    ##    ##     ##  ##     ##          ##   ##       ##       ##        
%  ######     ##    #### ##     ##         ##    ######## ##       ##        


l = plot(ax(2),nanmean(PID(25e3:45e3,:)),x_odor_LFP,'k+');
x = nanmean(PID(25e3:45e3,:));
y = x_odor_LFP; 
rm_this = isnan(x) | isnan(y);
x(rm_this) = []; y(rm_this) = [];
[rho,p] = corr(x(:),y(:),'type','Spearman');
legend(l,['\rho = ' oval(rho) ', p = ' oval(p)],'Location','southeast')

xlabel(ax(2),'Mean Stimulus (V)')
ylabel(ax(2),'Lag (ms)')
title(ax(2),'Odor \rightarrow LFP')

% ##       ######## ########     ##       ######## #### ########  #### ##    ##  ######   
% ##       ##       ##     ##     ##      ##        ##  ##     ##  ##  ###   ## ##    ##  
% ##       ##       ##     ##      ##     ##        ##  ##     ##  ##  ####  ## ##        
% ##       ######   ########        ##    ######    ##  ########   ##  ## ## ## ##   #### 
% ##       ##       ##             ##     ##        ##  ##   ##    ##  ##  #### ##    ##  
% ##       ##       ##            ##      ##        ##  ##    ##   ##  ##   ### ##    ##  
% ######## ##       ##           ##       ##       #### ##     ## #### ##    ##  ######   

l = plot(ax(3),nanmean(PID(25e3:45e3,:)),x_LFP_fA,'k+');
x = nanmean(PID(25e3:45e3,:));
y = x_LFP_fA; 
rm_this = isnan(x) | isnan(y);
x(rm_this) = []; y(rm_this) = [];
[rho,p] = corr(x(:),y(:),'type','Spearman');
legend(l,['\rho = ' oval(rho) ', p = ' oval(p)],'Location','southeast')

xlabel(ax(3),'Mean Stimulus (V)')
ylabel(ax(3),'Lag (ms)')
title(ax(3),'LFP \rightarrow Firing Rate')

%  ######  ######## #### ##     ##    ##       ######## #### ########  #### ##    ##  ######   
% ##    ##    ##     ##  ###   ###     ##      ##        ##  ##     ##  ##  ###   ## ##    ##  
% ##          ##     ##  #### ####      ##     ##        ##  ##     ##  ##  ####  ## ##        
%  ######     ##     ##  ## ### ##       ##    ######    ##  ########   ##  ## ## ## ##   #### 
%       ##    ##     ##  ##     ##      ##     ##        ##  ##   ##    ##  ##  #### ##    ##  
% ##    ##    ##     ##  ##     ##     ##      ##        ##  ##    ##   ##  ##   ### ##    ##  
%  ######     ##    #### ##     ##    ##       ##       #### ##     ## #### ##    ##  ######   

l = plot(ax(5),nanmean(PID(25e3:45e3,:)),x_odor_fA,'k+');
x = nanmean(PID(25e3:45e3,:));
y = x_odor_fA; 
rm_this = isnan(x) | isnan(y);
x(rm_this) = []; y(rm_this) = [];
[rho,p] = corr(x(:),y(:),'type','Spearman');
legend(l,['\rho = ' oval(rho) ', p = ' oval(p)],'Location','southeast')

xlabel(ax(5),'Mean Stimulus (V)')
ylabel(ax(5),'Lag (ms)')
title(ax(5),'Odor \rightarrow Firing Rate')

% ##       ####  ######   ##     ## ########          #######  ##    ## ##       ##    ## 
% ##        ##  ##    ##  ##     ##    ##            ##     ## ###   ## ##        ##  ##  
% ##        ##  ##        ##     ##    ##            ##     ## ####  ## ##         ####   
% ##        ##  ##   #### #########    ##    ####### ##     ## ## ## ## ##          ##    
% ##        ##  ##    ##  ##     ##    ##            ##     ## ##  #### ##          ##    
% ##        ##  ##    ##  ##     ##    ##            ##     ## ##   ### ##          ##    
% ######## ####  ######   ##     ##    ##             #######  ##    ## ########    ##    


p = '/local-data/DA-paper/Chrimson/scaled-light-flicker/v2';
[~, ~, fA, paradigm, orn, fly, AllControlParadigms] = consolidateData(p,true);

x = [.75:0.05:1.1 1.2:.1:3.6 3.8 4 4.2 4.5];
x = [0:0.1:0.7 x];
y = [0 4 12.8 24.1 36.2 48.1 60 70 95 114 134 151 167 184 201 219 236 252 269 283 299 314 328 343 357 370 383 394 404 413 419 422 422 421 419 418 417];
y = [0*(0:0.1:0.7) y ];

light_power_fit = fit(x(:),y(:),'spline');

% make new vectors for mean and range of the stimulus
a = 25e3; z = 55e3;
LED = NaN*fA;
LED_voltage = LED;
light_s = NaN*paradigm; light_m = NaN*paradigm;
for i = 1:length(paradigm)
	temp = AllControlParadigms(paradigm(i)).Name;
	LED_voltage(:,i) = (AllControlParadigms(paradigm(i)).Outputs(1,1:10:end));
	LED(:,i) = light_power_fit(AllControlParadigms(paradigm(i)).Outputs(1,1:10:end));
end

% remove some junk
rm_this = sum(fA) == 0;
LED(:,rm_this) = [];
fA(:,rm_this) = [];
paradigm(rm_this) = [];
orn(rm_this) = [];
fly(rm_this) = [];
LED_voltage(:,rm_this) = [];

x_LED_fA = NaN(1,width(LED));

for i = 1:width(LED)
	s = LED(25e3:45e3,i)-mean(LED(25e3:45e3,i)); s = s/std(s);
	f = fA(25e3:45e3,i)-mean(fA(25e3:45e3,i)); f = f/std(f);

	temp = xcorr(f,s);
	[~,x_LED_fA(i)] = max(temp(19.5e3+1:20.5e3));

end

x_LED_fA = x_LED_fA - 500;
x_LED_fA(x_LED_fA<-100) = NaN;

l = plot(ax(6),nanmean(LED(25e3:45e3,:)),x_LED_fA,'k+');
x = nanmean(LED(25e3:45e3,:));
y = x_LED_fA; 
rm_this = isnan(x) | isnan(y);
x(rm_this) = []; y(rm_this) = [];
[rho,p] = corr(x(:),y(:),'type','Spearman');
legend(l,['\rho = ' oval(rho) ', p = ' oval(p)],'Location','southeast')

xlabel(ax(6),'Mean Stimulus (\muW)')
ylabel(ax(6),'Lag (ms)')
title(ax(6),'Light \rightarrow Firing Rate')

set(ax(2),'YLim',[0 250],'XLim',[0 2])
set(ax(3),'YLim',[-100 0],'XLim',[0 2])
set(ax(5),'YLim',[0 100],'XLim',[0 2])
set(ax(6),'YLim',[0 100],'XLim',[0 200])

o = imread('../images/modules.png');
axes(ax(1));
imagesc(o);
axis ij
axis image
axis off

legend('boxoff')
prettyFig('fs=18;')

if being_published	
	snapnow	
	delete(gcf)
end


%% Version Info
%
pFooter;


