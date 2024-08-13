clear all
close all
clc

SNR = (0:0.2:3);
f = 0.5*erfc(sqrt(SNR));
subplot(2,1,1)
plot(SNR,log10(f))
grid on
subplot(2,1,2)
semilogy(SNR, f)
xlabel('SNR')
ylabel('Potência')
title('Plot com eixo Y em escala logarítmica')
grid on
