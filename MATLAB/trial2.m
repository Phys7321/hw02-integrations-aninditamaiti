clear all

dt = @(x,y) (4 ./sqrt(2)) ./sqrt(y .^4 - x .^4);
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