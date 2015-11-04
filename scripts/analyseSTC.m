% analyseSTC.m
% 
% created by Srinivas Gorur-Shandilya at 2:17 , 03 November 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [] = analyseSTC(S,R,K)

figure('outerposition',[0 0 1200 900],'PaperUnits','points','PaperSize',[1200 900]); hold on
subplot(3,1,1)
plot(1e-4*(1:length(S)),S)
title('Stimulus')

subplot(3,1,2)
raster2(R)
title('Response')

subplot(3,3,7), hold on
plot(K,'k')

% compute the STA
l(1) = plot(STA(R,S),'r');

% compute the STC
[K_STC,eigen_values] = STC(R,S);

l(2) = plot(K_STC(:,1),'b');
legend(l,{'STA','1st STC'})

subplot(3,3,8)
plot(K_STC(:,2),'b')
title('2nd STC')

subplot(3,3,9), hold on
plot(eigen_values,'k-+')
title('First six Eigenvalues')
