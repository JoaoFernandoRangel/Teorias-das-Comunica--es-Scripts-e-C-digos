

plot_Manchester = 0;

% Parâmetros
fs = 2e3; % Frequência de amostragem, define a frequência máxima a ser analisada no espectro de frequência
nb = 10; % Número de bits/símbolos
sb = randi([0 1], 1, nb) % Geração de sequência de bits aleatórios
% sb = [1 0 1 1 0 1 0 0];
nab = 8; % Número de amostras por bit/símbolo
A = 1; % Amplitude do sinal

% Geração do sinal Manchester
manchester = zeros(1, nb * nab); % Inicialização do vetor Manchester
for i = 1:nb
    if sb(i) == 0
        manchester((i-1)*nab + 1:i*nab) = [A * ones(1, nab/2), -A * ones(1, nab/2)];
    else
        manchester((i-1)*nab + 1:i*nab) = [-A * ones(1, nab/2), A * ones(1, nab/2)];
    end
end
subplot(2,1,1)
stem(sb);
grid on
subplot(2,1,2)
stem(manchester);
ylim([-1.2 1.2]);
grid on

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