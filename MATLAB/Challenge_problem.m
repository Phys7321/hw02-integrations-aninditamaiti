clear all

dt = @(x1,a) (2*sqrt(2)) ./sqrt(a .^4 - x1 .^4);
Amp = linspace(0.1,7,100);
T = [];
for i=1:length(Amp)
   Amp1 = Amp(i);
   x = linspace(0.01,Amp1*99/100);
   dT = trapz(x,dt(x,Amp1));
   T = [T, dT]
    
end

plot(Amp,T)
xlabel('Amplitude')
ylabel('Period')
title('Anharmonic oscillator')