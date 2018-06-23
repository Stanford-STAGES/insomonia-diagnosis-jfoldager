function freqRes = makefreqres(fRange, fBands, coeffs)
delta               = sigmoid(fRange,0.01,20);
tail                = fliplr(delta(delta<0.999));
a = find(fRange > fBands(2)); a = a(1);
delta(a:a + length(tail) - 1) = tail; delta(a + length(tail) - 1:end) = zeros(1,length(delta(a + length(tail) - 1:end)));
delta = coeffs(1)*delta;

theta               = sigmoid(fRange,fBands(2),10);
tail                = fliplr(theta(theta<0.999));
a = find(fRange > fBands(3)); a = a(1);
theta(a:a + length(tail) - 1) = tail; theta(a + length(tail) - 1:end) = zeros(1,length(theta(a + length(tail) - 1:end)));
theta = coeffs(2)*theta;

alpha               = sigmoid(fRange,fBands(3),10);
tail                = fliplr(alpha(alpha<0.999));
a = find(fRange > fBands(4)); a = a(1);
alpha(a:a + length(tail) - 1) = tail; alpha(a + length(tail) - 1:end) = zeros(1,length(alpha(a + length(tail) - 1:end)));
alpha = coeffs(3)*alpha;

beta               = sigmoid(fRange,fBands(4),10);
tail                = fliplr(beta(beta<0.999));
a = find(fRange > fBands(5)); a = a(1);
beta(a:a + length(tail) - 1) = tail; beta(a + length(tail) - 1:end) = zeros(1,length(beta(a + length(tail) - 1:end)));
beta = coeffs(4)*beta;

freqRes = [delta;theta;alpha;beta];
end