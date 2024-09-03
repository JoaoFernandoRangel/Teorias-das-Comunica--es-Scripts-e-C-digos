clear all
close all
clc
%{
Modulador e demodulador bpsk
%}
nb = 10;
seqb = randi([0 1], 1, nb);
% seqb = [0 0 0 1 0 0 1 0 1 0];
rb = 10e3; %Banda 
fc = 20e3; %Frequência da portadora
fs = 200e3; %frequencia de amostragem
nab = fs/rb; %amostras por bit
sinal = reshape(repmat(seqb, nab,1),1,nab*nb);
sinal_nrz = sinal*2-1;
ts = 1/fs;
t = 0:ts:(nab*nb-1)*ts;
phi = 2*pi*fc*t;
onda_transmitida = sinal_nrz.*cos(2*pi*fc*t);
subplot(1,2,1)
plot(t,onda_transmitida)
grid on
%Chegou o sinalrx
SNR = 10;
sinalrx = awgn(onda_transmitida, SNR); % + ruido gaussiano
snrzrx = sinalrx.*cos(2*pi*fc*t);
seqsrx = sum(reshape(snrzrx,nab,nb));
seqbrx = seqsrx>0;
subplot(1,2,2)
stem(seqbrx)
grid on