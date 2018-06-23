startup
L = 200; % number of data points
% generate data for example
C = [1 2; 1 1];
X = C * randn (2 , L );

figure, plot ( X (1 ,:) , X (2 ,:) , 'b.')
xlabel('x'),ylabel('y')
grid on, box off
axis ([ -1 1 -1 1]*8)

A = X * X';
[E , D ] = eig ( A );
err = A - E * D * E';
max ( abs ( err (:)));
d = diag ( D );
[ tmp , k ] = sort ( - d );
d = d ( k );
D = diag ( d );
E = E (: , k );
err = A - E * D * E';
max ( abs ( err (:)));
P = E';
Y = P * X ;
figure, plot ( Y (1 ,:) , Y (2 ,:) , 'b.')
axis ([ -1 1 -1 1]*8)
xlabel('PC 1'),ylabel('PC 2')
grid on, box off