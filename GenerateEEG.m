startup
c = [1,1,1,1];
N = 10000;
y = randn(1,N);
fs = 250;
T = 1/fs;
L = length(y);
t= linspace(0,L,L)*T;
% plot(t,y)
f = fs*linspace(0,floor(L/10),floor(L/10))/L;
Y = fft(y);
cutoff = [0.5,4.0,7.0,12.0,30.0];
powerCoefficients = [.1,.1,.1,.1];
filtered = [];
b = [];
for band = 1:length(cutoff)-1
 wl = 2*cutoff(band)/fs*pi;
 wh = 2*cutoff(band+1)/fs*pi;
 M = 512;
 bn = zeros(M);
 
 for i = 1:M    
  n = i-M/2;    
  if n == 0
   bn(i) = wh/pi - wl/pi;
  else
   bn(i) = (sin(wh*n))/(pi*n) - (sin(wl*n))/(pi*n);   
  end
 end
 bn = bn*kaiser(M,5.2);  
 
 b = [b bn];
 
 [w,h]=freqz(bn,1);
 filtered = [filtered conv(powerCoefficients(band)*1/(wh-wl)*bn, y)]; 
end
% plot(abs(fft(bn)))

fig = figure(1);
subplot(2, 1, 1)
plot(t,.1*y)
for i = 1:length(filtered)
  y_p = filtered(i);
  plot(t,y_p(floor(M/2):L+floor(M/2))+(i+1)*3)
end

y_t = y*0;
for i = 0:length(filtered)
    y_p = filtered[i]
    y_t = y_t + y_p[M//2:L+M//2]

plot(t,y_t-3)

axis('tight')
title('Time domain')
xlabel('Time (seconds)')
#plt.show()
    
#figure(figsize=[16, 10])
subplot(2, 1, 2)
plot(f,2*abs(Y[0:L//10]))
for i in range(0, len(filtered)):
  Y = filtered[i]
  Y = fft(Y [ M//2:L+M//2])
  plot(f,2*abs(np.real(Y[0:L//10])))
axis('tight')
legend(['synthetic','delta band, 0-4 Hz','theta band, 4-7 Hz','alpha band, 7-12 Hz','beta band, 12-30 Hz'])

for i in range(0, len(filtered)):   # plot filter's frequency response
  H = abs(fft(b[i], L))
  H = H*1.2*(max(Y)/max(H))
  plot(f, 3*np.real(H[0:L//10]), 'k')    
title('Frequency domain')
xlabel('Frequency (Hz)')
subplots_adjust(left=0.04, bottom=0.04, right=0.99, top=0.97)
axis('tight')





