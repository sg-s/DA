% stabiliseImageSequence.m
% stabilizes an image sequence using image registration
% 
% created by Srinivas Gorur-Shandilya at 1:32 , 02 December 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [stabilised_images] = stabiliseImageSequence(images)

reference_image = squeeze(images(:,:,1));

stabilised_images = images;
[o,m]=imregconfig('monomodal');

for i = 2:size(images,3)
	disp(i)
	b = squeeze(images(:,:,i));
	[tform] = imregtform(reference_image,b,'rigid',o,m);
	stabilised_images(:,:,i) = imwarp(b,tform,'OutputView',imref2d(size(b)));
end	