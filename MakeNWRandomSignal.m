% MakeNWRandomSignal
% this script makes a random binary signal according to what Nagel and Wilson said they did in their paper
% 
% created by Srinivas Gorur-Shandilya at 1:32 , 22 April 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

x = rand(1000,1);
x = repmat(x,1,50)';
x = x(:);
x(x<0.5)=0;
x(x>0)=1;

K = filter_exp(3e3,1,1:15e3); % 15 seconds long
K = K/K(1);

x = filter(K,1,x-.5);
x(x<0) = 0;
x(x>0)=1;