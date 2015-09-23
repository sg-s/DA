% MakeNWRandomSignal
% this script makes a random binary signal according to what Nagel and Wilson said they did in their paper
% 
% created by Srinivas Gorur-Shandilya at 1:32 , 22 April 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

% assume 1ms sampling time
x = rand(1000,1);
x = repmat(x,1,50)'; % this gives 20Hz
x = x(:);
x(x<0.5)=0;
x(x>0)=1;

% filter
K = filter_exp(3e3,1,1:15e3); % 15 seconds long filter, with tau=3s
K = K/K(1);
x = filter(K,1,x-.5);

% round
x(x<0) = 0;
x(x>0)=1;