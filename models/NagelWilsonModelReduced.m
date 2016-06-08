% NagelWilsonModelReduced.m
% specification of the Nagel-Wilson model for receptor binding, channel opening, and adptation via diffusable factor
% this is reduced from their formulation to include a conservation equation for the reaction species, which their model lacks
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function dy = NagelWilsonModelReduced(t,y,time,odor,p)

dy = NaN*y;

% calculate the odor at the time point
O = interp1(time,odor,t); % Interpolate the data set (ft,f) at time t

% receptor-ligand binding

R = y(1);
RX = y(2);
ORX = y(3);

dR = -R*(p.ka*p.sa + O*p.kb*p.sb) + RX*p.sa + p.sb*(1 - R - RX - ORX);
dRX = R*(p.ka*p.sa) - RX*(p.sa + O*p.theta*p.kb*p.sb) + ORX*p.sb;
dORX = RX*O*p.theta*p.kb*p.sb + (1 - R - RX - ORX)*p.theta*p.ka*p.sa - ORX*(p.sa+p.sb);

dy(1:3) = [dR dRX dORX];
  
% diffusible factor and channel opening
dy(4) = (RX+ORX)*p.ko*(1/y(5))*(1-y(4)) - p.kc*y(4); % according to Nagel and Wilson
dy(5) = p.A*y(4) - p.B*y(5);

