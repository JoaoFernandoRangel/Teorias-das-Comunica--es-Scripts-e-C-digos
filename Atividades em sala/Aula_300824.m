clear all
close all
clc
Nsim = 1e4; %Numero de simbolos
M = 16; %Símbolo por amostra
Ssimb = randi([0, M-1],1,Nsim);%Sequencia de simbolos
Sdmod_psk = pskmod(Ssimb,M); %Sequencia de dados modulada em psk
Sdmod_qam = qammod(Ssimb, M); %Sequencia de dados modulada em qam

% figure 1
% subplot 121
% plot(real(Sdmod_psk), imag(Sdmod_psk), '*r');
% ylim([-1.1 1.1])
% xlim([-1.1 1.1])
% grid on
% subplot 122
% plot(real(Sdmod_qam), imag(Sdmod_qam), '*b');
% ylim([-3.1 3.1])
% xlim([-3.1 3.1])
% grid on

Rs = 1e3;    %1000 bauds por segundo
fc = 2e3;    %Frequência da portadora
Fs = 4*fc;   %Frequência de amostragem
nas = Fs/Rs; %Numero de amostra por símbolo
Ts = 1/Fs; %Período entre amostras
Ns = Nsim*nas; %Número de amostras
df = Fs/Ns; %Incremento do vetor frequência
freq = -Fs/2:df:(Fs/2-df); %Vetor de frequência

seq_psk = reshape(repmat(Sdmod_psk,nas,1),1,Ns);
seq_qam = reshape(repmat(Sdmod_qam,nas,1),1,Ns);
sinal_psk = modulate(real(seq_psk), fc, Fs,'qam',imag(seq_psk));
sinal_qam = modulate(real(seq_qam), fc, Fs, 'qam', imag(seq_qam));
sinal_qam_prop = Qam_mod(real(seq_qam),fc,Fs,imag(seq_qam));
% Definir a taxa de amostragem
Fs = 10000; % Exemplo: 10 kHz
% Definir as frequências de corte
f1 = 1200; % 1.5 kHz
f2 = 2800; % 2.5 kHz
% Criar o filtro passa faixa
d = designfilt('bandpassiir', 'FilterOrder', 80, ...
               'HalfPowerFrequency1', f1, 'HalfPowerFrequency2', f2, ...
               'SampleRate', Fs);
% Exibir a resposta em frequência do filtro
fvtool(d, 'Fs', Fs);
y = filter(d,sinal_qam);
% figure 2
subplot 211
plot(freq,fftshift(abs(fft(y))));
grid on
subplot 212
plot(freq,fftshift(abs(fft(sinal_qam))));
grid on

