function [a,b,angle] = finding_length(data)
% Calculate the eigenvectors and eigenvalues
covar       = cov(data); [eigvec, eigval ] = eig(covar);
eigval      = diag(eigval);
% Get the largest eigenvector and eigenvalue
[maxEigVal, maxEigValIdx] = max(eigval);
maxEigVec   = eigvec(:, maxEigValIdx);
% Smallest eigenvector and eigenvalue
[minEigVal, minEigVecIdx] = min(eigval);
minEigVec   = eigvec(:,minEigVecIdx);
% Angle between the x-axis and the largest eigenvector
angle       = atan2(maxEigVec(2), maxEigVec(1));
% Shift it such that the angle is between 0 and 2pi
if(angle < 0),angle=angle+2*pi;if angle > pi, angle = angle - pi; end, end
mu = mean(data(:,1));
% Get 95% ci 
chisquare_val   = 2.4477;
a               = 2*chisquare_val*sqrt(maxEigVal);
b               = 2*chisquare_val*sqrt(minEigVal);
end