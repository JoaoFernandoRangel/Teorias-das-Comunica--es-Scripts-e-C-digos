clear all
close all
clc
A = 1;
Nb = 9e3 ; %Numero de bits
nab = 16; %numero de amostras de bit
Ns = nab*Nb; % 
Rb = 1e3; % Taxa de bits por segundo
Fs = nab*Rb; % Frequencia de amostragem
Ts = 1/Fs; %Período de amostragem
Ti = 0;
Tf = Nb/Rb; % Tempo final de transmissão
t = Ti:Ts:(Tf-Ts); %Vetor de tempo
df = Fs/Ns;
f = -(Fs/2):df:((Fs/2)-df);
sbits = randi([0 1], 1, Nb);
snrzu = A * reshape(repmat(sbits,nab,1),1,Ns);
% Como funciona o repmat : vetor, numero de vezes que vai repetir a linha, 1 vai definir quantas
% vezes a sequencia será repetida na linha
snrz = 2*snrzu - A;
fc1 = ones(1,nab);
fc2 = -ones(1,nab);
% Filtros casados para bit 1 e zero
SNR = -3; % Em dB
sinalrx = awgn(snrz,SNR);
% Sinal passa por um canal com ruído Gaussiano
sbrx = sinalrx(nab/2:nab:end);
lambda = 0;
sbitsrx = sbrx>lambda;
[bersf, rtsf] = biterr(sbits,sbitsrx);
% Recepção com filtro
sbcfrx  = sum(reshape(sinalrx,nab,Nb));
sbitscfrx = sbcfrx>lambda;
[bercf, rtcf] = biterr(sbits,sbitscfrx);
%{
bersf
rtsf
bercf
rtcf
%}

%Gerando M-PAM
M = 8;
n = log2(M);
sbq = randi([0,M-1],1,Nb/n);
seqsb = reshape(repmat(sbq,nab*n,1),1,Nb*nab);
s4pam = pammod(seqsb,M);
subplot(2,1,1)
plot(real(s4pam), imag(s4pam), '.r')
grid on
subplot(2,1,2)
plot(real(s4pam))
xlim([0 400]);
grid on
close all
plot(fftshift(abs(fft(s4pam))))



