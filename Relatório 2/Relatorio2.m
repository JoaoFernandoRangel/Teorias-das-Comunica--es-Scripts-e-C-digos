clear all
close all
clc
%Variáveis de controle do código
plotaGraficos = 0;


% Parâmetros
nb = 100; % Número de bits/símbolos
sb = randi([0 1], 1, nb); % Geração de sequência de bits aleatórios
nab = 8; % Número de amostras por bit/símbolo
A = 1; % Amplitude do sinal

% NRZ Unipolar
sNRZ_unipolar = A * (reshape(repmat(sb, nab, 1), 1, nab * nb));

% NRZ Bipolar
sb_bipolar = 2*sb - 1; % Mapeia 0 -> -1, 1 -> 1
sNRZ_bipolar = A * (reshape(repmat(sb_bipolar, nab, 1), 1, nab * nb));

% RZ Unipolar
bit1 = [1 1 1 1 0 0 0 0]; % Padrão para RZ
sRZ_unipolar = sNRZ_unipolar .* repmat(bit1, 1, nb);

% RZ com Inversão do 1
sRZ_inversao = sNRZ_bipolar .* repmat(bit1, 1, nb);

% Manchester
sManchester = zeros(1, 2 * nab * nb);
for i = 1:nb
    if sb(i) == 0
        sManchester((2*i-1)*nab-nab+1:(2*i-1)*nab) = A * ones(1, nab);
        sManchester((2*i-1)*nab+1:2*i*nab) = -A * ones(1, nab);
    else
        sManchester((2*i-1)*nab-nab+1:(2*i-1)*nab) = -A * ones(1, nab);
        sManchester((2*i-1)*nab+1:2*i*nab) = A * ones(1, nab);
    end
end

% 4-PAM
M0 = 4; % Número de níveis (4-PAM)
levels = [-3, -1, 1, 3]; % Definição dos níveis 4-PAM

% Convertendo os bits para símbolos 4-PAM
sb_reshaped = reshape(sb, 2, []); % Agrupando os bits em pares
symbols = sb_reshaped(1,:) * 2 + sb_reshaped(2,:) + 1; % Convertendo bits para índice de nível
sPAM4 = levels(symbols); % Mapeando para os níveis 4-PAM

% Geração do sinal 4-PAM
sPAM4_signal = reshape(repmat(sPAM4, nab, 1), 1, nab * length(sPAM4));

% 8-PAM
sb_pam8 = randi([0 7], 1, nb); % Geração de sequência de símbolos 8-PAM
M = 8; % Número de níveis (8-PAM)
pam8_levels = [-7, -5, -3, -1, 1, 3, 5, 7]; % Definição dos níveis 8-PAM
sPAM8 = pam8_levels(sb_pam8 + 1); % Mapeamento dos símbolos para níveis 8-PAM
sPAM8_signal = reshape(repmat(sPAM8, nab, 1), 1, nab * nb);

%%
clear all
close all
clc
fs = 40; % Frequência de amostragem
T = 1/fs; % Período de amostragem
N = 40; % Número de amostras
t = (0:N-1)*T'; % Vetor de tempo
xn = cos(2*pi*3*t)+cos(2*pi*15*t); % Sinal de entrada

% Plot do sinal
subplot(2,1,1)
plot(t, xn)
title('Sinal no Tempo')
xlabel('Tempo (s)')
ylabel('Amplitude')
grid on
% FFT
y = fft(xn);
f = (0:N-1)*(fs/N); % Eixo de frequência
% Plot do espectro
subplot(2,1,2);
stem(f, abs(y))
title('Espectro de Frequência')
xlabel('Frequência (Hz)')
ylabel('Magnitude')
xlim([0 fs/2]) % Limita o eixo x para melhor visualização
grid on

%%

if(plotaGraficos == 1)
% Plotagem
figure (1);
subplot(2,1,1);
plot(sNRZ_unipolar);
xlim([0 200]);
ylim([-0.1 1.2]);
grid on
title('NRZ Unipolar');
xlabel('Tempo');
ylabel('Amplitude');
subplot(2,1,2);
plot(sNRZ_bipolar);
xlim([0 200]);
ylim([-1.2 1.2]);
grid on
title('NRZ Bipolar');
xlabel('Tempo');
ylabel('Amplitude');
figure (2);
subplot(2,1,1);
plot(sRZ_unipolar);
xlim([0 200]);
ylim([-0.1 1.2]);
grid on
title('RZ Unipolar');
xlabel('Tempo');
ylabel('Amplitude');
subplot(2,1,2);
plot(sRZ_inversao);
xlim([0 200]);
ylim([-1.2 1.2]);
grid on
title('RZ com Inversão do 1');
xlabel('Tempo');
ylabel('Amplitude');
figure (3);
plot(sManchester);
xlim([0 200]);
ylim([-1.2 1.2]);
grid on
title('Manchester');
xlabel('Tempo');
ylabel('Amplitude');
figure (4)
subplot(2,1,1)
plot(sPAM4_signal);
xlim([0 200]);  
ylim([-4.2 4.2]);
grid on 
title('4-PAM');
xlabel('Tempo');
ylabel('Amplitude');
subplot(2,1,2)
plot(sPAM8_signal);
xlim([0 200]);
grid on
title('8-PAM');
xlabel('Tempo');
ylabel('Amplitude');
end
