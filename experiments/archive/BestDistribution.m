% BestDistribution.m
% finds best control signal distribution to achieve target stimulus distribution 
% 
% assumptions:
% 
% - global variable p1 corresponding to a structure for pDeliverySystem (control -> MFC)
% - global variable p2 corresponding to a structure for pDeliverySystem (MFC -> PID)
% - data handed to FitModel2Data should have the following fields:
% - data.stimulus: (blank)
% - data.response: the target distribution we want to achieve, whose domain is px (see below)
%
% BestDistribution is meant to be optimised by FitModel2Data. To understand how FitModel2Data works, see its documentation first. BestDistribution does not require an input stimulus. 
% 
% Every time BestDistribution is called (by FitModel2Data), it does the following:
% 
% - using input argument p, generate a control signal distribution over domain cx using dist_gamma2
% - sample deterministically from this using pdfrnd and frozen noise
% - feed this control to pDeliverySystem using p1 and get the output
% - feed that output as input to pDeliverySystem, using p2 and getting the output
% - histogram it on domain px
% - return this histogram as py so that FitModel2Data can penalise appropriately 
% 
% created by Srinivas Gorur-Shandilya at 1:53 , 26 February 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [py,MFC_pred,PID_pred] = BestDistribution(~,p)

global p1
global p2


% some critical parameters
T = 20; 		% how long is your stimulus
dt = 1e-3; 		% sampling time for this simulation
tc = .05; 		% switching time
cx = 0:1e-2:5;  % control signal domain
px = 0:1e-2:5;  % PID domain
MFC_min = 5/200;% 200:1 turndown ratio

% these bounds are the same as in dist_gamma2, which is the core of this function
lb.A = 0;
lb.n1 = 0;
lb.n2 = 0;
lb.xmin = 0;
lb.xmax = 0;
ub.xmin = 1;
ub.xmax = 1;
ub.n1 = 10;
ub.n2 = 10;
lb.x1 = 0;
lb.x2 = 0;
ub.x1 = 5;
ub.x2 = 5;

% make parameters readable to FitModel2Data
p.xmin;
p.xmax;
p.n1;
p.n2;
p.A;
p.x1;
p.x2;




% generate distribution
cy = dist_gamma2(cx,p);
cy(isnan(cy)) = 0;

% sample from distribution 
n = floor(T/tc);
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',1984)); 
s = repmat(pdfrnd(cx,cy,n),1,floor(tc/dt));
s = s';
s = s(:);
s(s<MFC_min) = MFC_min;

py=0*px;

% get pid prediction
MFC_pred = pDeliverySystem(s,p1);
MFC_pred = (100*MFC_pred)./(100*MFC_pred+2000);
PID_pred = pDeliverySystem(MFC_pred,p2);

% get pid distribution 
py = hist(PID_pred(1e3:end),px);
py(isnan(py)) = 0;
py = py/max(py);

% penalise bullshit distributions
py(py<eps) = 0;
if length(nonzeros(py)) < 10
    py = 1*py;
end

