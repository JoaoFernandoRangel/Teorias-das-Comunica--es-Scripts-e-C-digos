%%
%Exemplo 4.1 - Filtro casado
nAmostras = 16;
sinal = [ones(1,8) (-1*ones(1,8))];
subplot(2,2,1)
plot(sinal);
ylim([-1.2 1.2]);
grid on
fc = -sinal;
subplot(2,2,2)
plot(fc);
ylim([-1.2 1.2]);
grid on
subplot(2,2,3);
srec = conv(sinal,fc);
plot(srec);
grid on
%%
%Exemplo 4.2
clear all
close all
clc


T = 16e-3;
Na = 8;
t = linspace(0,32e-3,15);
s1 = [1 1 1 1 -1 -1 -1 -1];
h1 = [-1 -1 -1 -1 1 1 1 1];
h2 = [-1 -1 1 1 1 1 -1 -1];
s2 = h2;
y11 = conv(s1,h1);
y12 = conv(s1,h2);
y21 = conv(s2,h1);
y22 = conv(s2,h2);
subplot(2,2,1)
plot(t, y11)
grid on
subplot(2,2,2)
plot(t, y12)
grid on
subplot(2,2,3)
plot(t, y21)
grid on
subplot(2,2,4)
plot(t, y22)
grid on
%%
%Exemplo real
clear all
close all
clc
nsb = 16;
numB = 1e2;
seqb = randi([0,1],1,numB);
A = 1;
SnrzB = A*(2*(reshape(repmat(seqb,nsb,1),1,nsb*numB))- 1);
SNR = 0;
sinalRecebido = awgn(SnrzB, SNR);
lambda = 0;
seqbr= sum(reshape(sinalRecebido,nsb,numB)) > lambda;
% subplot(2,1,1)
% stem(sequenciaBinaria)
% ylim([-0.1 1.1]);
% xlim([0 25]);
% grid on
% subplot(2,1,2)
% stem(seqbr)
% ylim([-0.1 1.1]);
% xlim([0 25]);
% grid on
seqbt = logical(seqb);
biterr(seqbt, seqbr)
bit1 = repmat([ones(1,nsb/2) -1*ones(1,nsb/2)],1,numB);
bit0= repmat([-1*ones(1,nsb/4) ones(1,nsb/2) -1*ones(1,nsb/4)],1,numB);
fc1 = -reshape(bit1,nsb,numB);
fc0 = reshape(bit0,nsb,numB);
snrzu1 = (reshape(repmat(seqb, nsb, 1),1,nsb*numB));
snrzu0 = -snrzu1+1;
sinal = snrzu1.*bit1 + snrzu0.*bit0;
% plot(sinal)
sinalr = awgn(sinal, SNR);
srconv = reshape(sinalr,nsb,numB);
sf1 = sum(conv(srconv, fc1));
sf0 = sum(conv(srconv, fc0));

















