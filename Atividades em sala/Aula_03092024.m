clear all
close all
clc
nb = 10;
seqb = randi([0 1], 1, nb);
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
plot(t,onda_transmitida)