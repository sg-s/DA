
%% The Nagel-Wilson Model
% In this document we try to make sense of the Nagel Wilson Model and see if it explains gain scaling and LFP slowdown.

pHeader;

%% 1: Using the Nagel-Wilson Model as is
% In this section, we use the Nagel-Wilson model as described in the paper. This has 6 ODEs, and the output of the model is the channel open probability. In the following figure, I plot the response of this model to a test pulse of stimulus. 

clear p
p.theta = 100;
p.   ka = 0.001;
p.   kb = 1;
p.   sa = 10;
p.   sb = 100;
p.   k0 = 500;
p.    A = 0.8000;
p.    B = 0.6000;
p.   ko = 500;
p.   kc = 10;

S = .1+zeros(10e3,1);
S(4e3:5e3) = 1;
[C,D,R] = NagelWilsonIntegrate(S,p);

figure('outerposition',[0 0 1200 700],'PaperUnits','points','PaperSize',[1200 700]); hold on
subplot(2,2,1), hold on
plot(S)
ylabel('Stimulus')
subplot(2,2,2), hold on
plot(R)
legend({'R','R*','OR','OR*'})
ylabel('Receptor Species')
subplot(2,2,3), hold on
plot(C)
ylabel('C_{open}')
subplot(2,2,4), hold on
plot(D)
ylabel('Diffusible factor')

prettyFig();

if being_published
  snapnow
  delete(gcf)
end

%%
% We see that there is a fatal flaw in the model: there is no conservation equation for the four reaction species. Thus, they are allowed to grow without bound, which is obviously wrong. 

%% 2. Nagel Wilson Model with conservation equation
% In this section, we add a conservation equation to the model:
% 
% $$ OR=R_{total}-R-R*-OR*$$
%

%%
% Doing so eliminates one of the ODEs. The ODEs now are:
%
% $$ \dot{R}=-R(k_{a}s_{a}+Ok_{b}s_{b})+R*s_{a}+(R_{total}-R-R*-OR*)s_{b} $$
%

%%
% and so on.

%%
% In the following figure, we plot the responses to a pulse of odorant.

S = .1+zeros(10e3,1);
S(4e3:5e3) = 1;
[C,D,R] = NagelWilsonIntegrate2(S,p);

figure('outerposition',[0 0 1200 700],'PaperUnits','points','PaperSize',[1200 700]); hold on
subplot(2,2,1), hold on
plot(S)
ylabel('Stimulus')
subplot(2,2,2), hold on
plot(R)
legend({'R','R*','OR','OR*'})
ylabel('Receptor Species')
subplot(2,2,3), hold on
plot(C)
ylabel('C_{open}')
subplot(2,2,4), hold on
plot(D)
ylabel('Diffusible factor')

prettyFig();

if being_published
  snapnow
  delete(gcf)
end

%%
% OK, this is somewhat reasonable. Now, we check to see if this model predicts a slowdown in response with increased stimulus. We do this by stimulating the model with pulses on top of increasing backgrounds. 

clear ax
figure('outerposition',[0 0 1200 700],'PaperUnits','points','PaperSize',[1200 700]); hold on
for i = 1:4
  ax(i) = subplot(2,2,i); hold on
end

background_levels = logspace(-2,0,5);
c = parula(length(background_levels)+1);
foreground_level = 2;

for i = 1:length(background_levels)
  S = background_levels(i) + zeros(10e3,1);
  S(5e3:6e3) = foreground_level;
  plot(ax(1),S,'Color',c(i,:));

  [C,D,R] = NagelWilsonIntegrate2(S,p);
  plot(ax(2),R(:,2)+R(:,4),'Color',c(i,:))

  plot(ax(3),C,'Color',c(i,:))

  plot(ax(4),D,'Color',c(i,:))
end

ylabel(ax(1),'Stimulus')
ylabel(ax(2),'R* + OR*')
ylabel(ax(3),'C_{open}')
ylabel(ax(4),'Diffusible factor')

prettyFig();

if being_published
  snapnow
  delete(gcf)
end


%%
% That looks interesting. What if we rescale the channel open probabilities to look more closely at the kinetics? 


clear ax
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
for i = 1:2
  ax(i) = subplot(1,2,i); hold on
end

for i = 1:length(background_levels)
  S = background_levels(i) + zeros(10e3,1);
  S(5e3:6e3) = foreground_level;
  plot(ax(1),S,'Color',c(i,:));

  [C,D,R] = NagelWilsonIntegrate2(S,p);

  C = C - mean(C(4e3:5e3));
  C = C/max(C(5e3:6e3));

  plot(ax(2),C,'Color',c(i,:))
end

ylabel(ax(1),'Stimulus')
ylabel(ax(2),'C_{open} (rescaled)')
set(ax(1),'XLim',[4000 8e3])
set(ax(2),'XLim',[4000 8e3],'YLim',[0 1])


prettyFig();

if being_published
  snapnow
  delete(gcf)
end


%%
% We clearly see that it slows down with increasing background stimulus. 


%% Version Info
%
pFooter;
