clear all
close all
clc
%Construindo o sinal
Ns = 10e3;
M = 8;
n = log2(M); %num de bits por simbolo
ssimb = randi([0, M-1],1,Ns);
Rs = 1e3; %taxa de símbolo
Ts = 1/Rs; %Tempo de símbolo
Rb = Rs*n; %Taxa de bits
Tb = 1/Rb; %Tempo de bit
fc = 5e3; % frequencia da portadora, tem que ser um multiplo inteiro de Rs
Fa = 25e3;%Frequência de amostragem
Ta = 1/Fa;%Período de amostragem
nas = Fa/Rs; %numero de amostra por simbolo
t = 0:Ta:Ta*(nas*Ns-1);
%Sinal construido
seqs = reshape(repmat(ssimb,nas,1),1,nas*Ns); % Sequencia de simbolos
%salada de reshape e repmat
seqsa = 2*seqs*pi./M; %theta de cada símbolo
Amp = 1; %amplitude do sinal
sinaltx = Amp*sin(2*pi*fc.*t + seqsa); %Sinal a ser transmitido
%Agora passamos o sinal por um AWGN
SNR = 3;
sinalrx = awgn(sinaltx, SNR, 'measured');
seqI = sinalrx .* cos(2*pi*fc.*t); %Sinal em fase
seqQ = sinalrx .* sin(2*pi*fc.*t); %Sinal em quadratura
seqarx = atan(seqQ./seqI); %Sequencia de angulo do Rx
spos = seqarx>0;
seqsarx = seqarx +(1-spos)*2*pi;
seqsrx = M*(seqsarx)./(2*pi);
ssimbrx = round(mean(reshape(seqsrx,nas,Ns))./nas); % Não sabe se ta certo
[ber, rt] = biterr(ssimb, ssimbrx)



