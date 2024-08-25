clear all
close all
clc

%Variáveis de controle do código
plotNRZ_Unipolar = 0;
plotNRZ_Bipolar =  0;
plotRZ_Unipolar =  0;
plotRZ_Inverso =   0;
plot_Manchester =  0;
plot_4_PAM =       0;
plot_8_PAM =       1;

% Parâmetros
fs = 2e3; % Frequência de amostragem, define a frequência máxima a ser analisada no espectro de frequência
nb = 10000; % Número de bits/símbolos
sb = randi([0 1], 1, nb); % Geração de sequência de bits aleatórios
nab = 8; % Número de amostras por bit/símbolo
A = 1; % Amplitude do sinal

% NRZ Unipolar
sNRZ_unipolar = A * (reshape(repmat(sb, nab, 1), 1, nab * nb));

% Plotagem do sinal e espectro
if (plotNRZ_Unipolar == 1)
    figure (1)
    % Plot do sinal NRZ Unipolar
    subplot(2,1,1);
    plot(sNRZ_unipolar);
    grid on;
    ylim([-0.1 1.2]);
    title('Sinal NRZ Unipolar');
    
    % Cálculo da Transformada de Fourier e espectro de potência
    N0 = length(sNRZ_unipolar);
    f0 = (0:N0-1)*(1/N0); % Frequências
    S0 = abs(fft(sNRZ_unipolar)).^2 / N0; % Espectro de potência
    
    % Plot do espectro de potência
    subplot(2,1,2);
    plot(f0*fs, S0);
    xlim([0 max(f0*fs)/2]);
    %     ylim([0 20]);
    grid on;
    title('Espectro de Potência do Sinal NRZ Unipolar');
    xlabel('Frequência (Hz)');
    ylabel('Densidade Espectral de Potência');
end

% NRZ Bipolar
sb_bipolar = 2*sb - 1; % Mapeia 0 -> -1, 1 -> 1
sNRZ_bipolar = A * (reshape(repmat(sb_bipolar, nab, 1), 1, nab * nb));
if (plotNRZ_Bipolar == 1)
    figure (2)
    % Plot do sinal NRZ Bipolar
    subplot(2,1,1);
    plot(sNRZ_bipolar);
    grid on;
    ylim([-1.2 1.2]);
    title('Sinal NRZ Bipolar');
    
    % Cálculo da Transformada de Fourier e espectro de potência
    N1 = length(sNRZ_bipolar);
    f1 = (0:N1-1)*(1/N1); % Frequências
    S1 = abs(fft(sNRZ_bipolar)).^2 / N1; % Espectro de potência
    
    % Plot do espectro de potência
    subplot(2,1,2);
    plot(f1*fs, S1);
    xlim([0 max(f1*fs)/2]);
    %     ylim([0 25]);
    grid on;
    title('Espectro de Potência do Sinal NRZ Bipolar');
    xlabel('Frequência (Hz)');
    ylabel('Densidade Espectral de Potência');
end

% RZ Unipolar
bit1 = [1 1 1 1 0 0 0 0]; % Padrão para RZ
sRZ_unipolar = sNRZ_unipolar .* repmat(bit1, 1, nb);
if (plotRZ_Unipolar == 1)
    figure (3)
    % Plot do sinal RZ Unipolar
    subplot(2,1,1);
    plot(sRZ_unipolar);
    grid on;
    ylim([-0.2 1.2]);
    title('Sinal RZ Unipolar');
    
    % Cálculo da Transformada de Fourier e espectro de potência
    N2 = length(sRZ_unipolar);
    f2 = (0:N2-1)*(1/N2); % Frequências
    S2 = abs(fft(sRZ_unipolar)).^2 / N2; % Espectro de potência
    
    % Plot do espectro de potência
    subplot(2,1,2);
    plot(f2*fs, S2);
    xlim([0 max(f2*fs)/2]);
    %     ylim([0 25]);
    grid on;
    title('Espectro de Potência do Sinal RZ Unipolar');
    xlabel('Frequência (Hz)');
    ylabel('Densidade Espectral de Potência');
end

% RZ com Inversão do 1
sRZ_inversao = sNRZ_bipolar .* repmat(bit1, 1, nb);
if (plotRZ_Inverso == 1)
    figure (4)
    % Plot do sinal RZ com inversão
    subplot(2,1,1);
    plot(sRZ_inversao);
    grid on;
    ylim([-1.2 1.2]);
    title('Sinal RZ com inversão');
    
    % Cálculo da Transformada de Fourier e espectro de potência
    N3 = length(sRZ_inversao);
    f3 = (0:N3-1)*(1/N3); % Frequências
    S3 = abs(fft(sRZ_inversao)).^2 / N3 ; % Espectro de potência
    
    % Plot do espectro de potência
    subplot(2,1,2);
    plot(f3*fs, S3);
    xlim([0 max(f3*fs)/2]);
    %     ylim([0 25]);
    grid on;
    title('Espectro de Potência do Sinal RZ com Inversão');
    xlabel('Frequência (Hz)');
    ylabel('Densidade Espectral de Potência');
end

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
if (plot_Manchester == 1)
    figure (5)
    % Plot do sinal RZ com inversão
    subplot(2,1,1);
    plot(sManchester);
    grid on;
    ylim([-1.2 1.2]);
    title('Sinal Manchester');
    
    % Cálculo da Transformada de Fourier e espectro de potência
    N4 = length(sManchester);
    f4 = (0:N4-1)*(1/N4); % Frequências
    S4 = abs(fft(sManchester)).^2 / N4 ; % Espectro de potência
    
    % Plot do espectro de potência
    subplot(2,1,2);
    plot(f4*fs, S4);
    xlim([0 max(f4*fs)/2]);
    grid on;
    title('Espectro de Potência do Sinal Manchester');
    xlabel('Frequência (Hz)');
    ylabel('Densidade Espectral de Potência');
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

if (plot_4_PAM == 1)
    figure (6)
    % Plot do sinal RZ com inversão
    subplot(2,1,1);
    plot(sPAM4_signal);
    grid on;
    title('Sinal 4-PAM');
    
    % Cálculo da Transformada de Fourier e espectro de potência
    N5 = length(sPAM4_signal);
    f5 = (0:N5-1)*(1/N5); % Frequências
    S5 = abs(fft(sPAM4_signal)).^2 / N5 ; % Espectro de potência
    
    % Plot do espectro de potência
    subplot(2,1,2);
    plot(f5*fs, S5);
    xlim([0 max(f5*fs)/2]);
    grid on;
    title('Espectro de Potência do 4-PAM');
    xlabel('Frequência (Hz)');
    ylabel('Densidade Espectral de Potência');
end

length

%8-PAM
sb_pam8 = randi([0 7], 1, nb) % Geração de sequência de símbolos 8-PAM
M = 8; % Número de níveis (8-PAM)
pam8_levels = [-7, -5, -3, -1, 1, 3, 5, 7]; % Definição dos níveis 8-PAM
sPAM8 = pam8_levels(sb_pam8 + 1); % Mapeamento dos símbolos para níveis 8-PAM
sPAM8_signal = reshape(repmat(sPAM8, nab, 1), 1, nab * nb);

if (plot_8_PAM == 1)
    figure (7)
    % Plot do sinal RZ com inversão
    subplot(2,1,1);
    plot(sPAM8_signal);
    grid on;
    title('Sinal 8-PAM');
    
    % Cálculo da Transformada de Fourier e espectro de potência
    N6 = length(sPAM8_signal);
    f6 = (0:N6-1)*(1/N6); % Frequências
    S6 = abs(fft(sPAM8_signal)).^2 / N6 ; % Espectro de potência
    
    % Plot do espectro de potência
    subplot(2,1,2);
    plot(f6*fs, S6);
    xlim([0 max(f6*fs)/2]);
    grid on;
    title('Espectro de Potência do 8-PAM');
    xlabel('Frequência (Hz)');
    ylabel('Densidade Espectral de Potência');
    
close all
end

