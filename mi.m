% mi.m
% 
% created by Srinivas Gorur-Shandilya at 11:13 , 21 August 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function m = mi(x,y)

nbins = 50;

assert(isvector(x) & isvector(y),'Inputs should be vectors')
assert(length(x) == length(y),'Inputs should have equal lengths')
assert(~any(isnan(x)) & ~any(isnan(y)),'Inputs should not have NaN')

x = x(:);
y = y(:);

% get both x and y from [0 1]
x = x -min(x); x = x/max(x);
y = y - min(y); y = y/max(y);

% compute all probability distributions
px = hist(x,nbins)/length(x);
py = hist(y,nbins)/length(x);

pxy = hist3([x y],[nbins nbins])/length(x);


m = 0;
for i = 1:nbins
	for j = 1:nbins
		temp = pxy(i,j)*log2(pxy(i,j)/(px(i)*py(j)));
		if ~isnan(temp)
			m = m+temp;
		end
	end
end
