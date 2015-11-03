function [V,D] = temp(S,spikes)
%Simulate the output of quad pair of Gabors
% rg   = (-7:8)';
% sigma = 1.6;
% gab1 = exp(-rg.^2/1/sigma^2).*cos(rg*2*pi/4);
% gab2 = exp(-rg.^2/1/sigma^2).*sin(rg*2*pi/4);
% gab1 = gab1 - mean(gab1)*exp(-rg.^2/1/sigma^2)/mean(exp(-rg.^2/1/sigma^2));
 
% %Probe filters with random stimuli
% X = randn(2000,16);
% %response w/Poisson noise
% r = poissrnd(.5*((X*gab1).^2+(X*gab2).^2));

r = spikes(:);
X = S;
 
%Recover the filter
%Create the matrix of products
X2 = anovize(X);
 
%Estimate the params through linear regression
ws = [X2,ones(size(X2,1),1)]\r;
 
%Form the Q matrix
Q = deanovize(ws,size(X,2));
 
%Take the first 6 eigenvalues
[V,D] = eigs(Q,6);
 
 
function [Y] = anovize(X)
	Y = zeros(size(X,1),(size(X,2)+1)*(size(X,2)/2));
	Y(:,1:size(X,2)) = X.^2;
	offset = size(X,2);
	for ii = 2:size(X,2)
		Y(:,offset + (1:size(X,2)+1-ii)) = 2*X(:,1:end-ii+1).*X(:,ii:end);
		offset = offset + size(X,2) + 1 - ii;
	end
end
 
function [Y] = deanovize(w,n)
	Y = zeros(n,n);
	offset = 0;
	for ii = 1:n
		rg = offset + (1:n+1-ii);
		Y(sub2ind([n,n],(1:n+1-ii)',(ii:n)')) = w(rg);
		offset = rg(end);
	end
	Y = Y + triu(Y)';

end

end