% freqResp.m
% computes the frequency response from raw, time series data
% 
% created by Srinivas Gorur-Shandilya at 1:31 , 10 June 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [H] = freqResp(x,y)

% compute the cross correlation function
xy = xcorr(x,y,'unbiased');

% compute the autocorrelation function of the input
xx = xcorr(x,x,'unbiased');

% compute the fourier transforms 
Sxy = fft(xy);
Sxx = fft(xx);

H = Sxy./Sxx;