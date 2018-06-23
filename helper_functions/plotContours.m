function plotContours(data,contours,options)

mu = [1 2;-3 -5];
sigma = cat(3,[2 0;0 .5],[1 0;0 1]);

mu = contours.mu(:,1:2); 
sigma = cat(3,contours.Sigma(1,1:2,:)); 
gm = gmdistribution(mu,sigma);
gmPDF = @(x,y)pdf(gm,[x y]);
ezcontour(gmPDF,[-10 10],[-10 10]);
title('Contour lines of pdf');




    figure,
    hold on
    h1 = scatter(data(:,1),data(:,2)); h = gca;
    mu = contours.mu; 
    sigma = cat(3,contours.Sigma); 
    p = ones(1,size(mu,1))/size(mu,1);
    gm = gmdistribution(mu,sigma,p);
    gmPDF = @(x,y)pdf(gm,[x y]);
    h2 = ezcontour(gmPDF,[min(data(:,1)) max(data(:,1))],[min(data(:,2)) max(data(:,2))],200);
    xlabel(options.xl)
    ylabel(options.yl)
    title(options.tl)
    grid on
    box on
    hold off
end