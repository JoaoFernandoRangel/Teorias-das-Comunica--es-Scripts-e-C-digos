clear all
close all
clc

N = 100; %Número de simbolos
M = 8; %M-PSK ou M-QQAM
nas = 8; %Amostras por simbolo em banda base
Rs = 1000; %Taxa de simbolos
fc = 4000; %Frequência da portadora
seqsimb = randi([0,M-1],1,N);
sibb = reshape(repmat(seqsimb,nas,1),1,nas*N); %Símbolos em banda base
angpsk = (2*pi*sibb)./M; %Angulos em PSK
simpsk = cos(angpsk) + 1i*sin(angpsk); %Simbolos em PSK
simqam = qammod(sibb,M); %Modulação em QAM

% A partir daqui é o modulador psk + qam
% close all 
clc
Fs = 4*(fc + Rs/2); %Frequência de amostragem
[P,Q] = rat(Fs/Rs); %Usa o numero inteiro mais próximo para a quantidade de amostrar na banda passante
sqam = resample(simqam,P,Q,10); 
spsk = resample(simpsk,P,Q,10);
Ts = 1/Fs;
na = length(sqam);
t = 0:Ts:Ts*(na-1);
sinalqamtx = real(sqam).*cos(2*pi*fc.*t) + imag(sqam).*sin(2*pi*fc.*t);
sinalpsktx = real(spsk).*cos(2*pi*fc.*t) + imag(spsk).*sin(2*pi*fc.*t);

%Canal de transmissão
SNR = -10:1:10;
for k = 1:length(SNR)
    sinalqamrx = awgn(sinalqamtx, SNR(k));
    sinalpskrx = awgn(sinalpsktx, SNR(k));   
    demodpsk_qam; %Demodulador psk e qam
    seqspsk = mean(reshape(spskrx,nas,N));
    angrx = atan(imag(seqspsk)./real(seqspsk));
    angrx(angrx<0) = 2*pi + angrx(angrx<0); % Ajuste de ângulo, soma uma volta no círculo
    simbpskrx = angrx*M./(2*pi);
    [bpsk(k),rtpsk(k)] = biterr(round(seqsimb), round(simbpskrx));
    seqsqamrx = mean(reshape(sqamrx,nas,N));
    simbqamrx = qamdemod(seqsqamrx, M);
    [bqam(k),rtqam(k)] = biterr(round(seqsimb), round(simbqamrx));
end

subplot(2,2,1)
plot(real(simpsk), imag(simpsk),'*b');
grid on
hold on
plot(real(simqam), imag(simqam),'*r');
subplot(2,2,2)
plot(real(seqspsk), imag(seqspsk),'*r');
hold on
plot(real(simpsk), imag(simpsk),'*b');
grid on
subplot(2,2,3)
plot(real(seqsqamrx), imag(seqsqamrx),'*r');
hold on
plot(real(simqam), imag(simqam),'*b');
grid on
subplot(2,2,4)
plot(SNR, rtpsk);
grid on










