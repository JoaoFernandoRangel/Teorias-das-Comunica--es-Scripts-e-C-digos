% Especificações do filtro
Fs = 100e6; % Frequência de amostragem (100 MHz, ajuste conforme necessário)
Fpass1 = 10e6; % Frequência de passagem baixa (10 MHz)
Fpass2 = 20e6; % Frequência de passagem alta (20 MHz)
Fstop1 = 5e6;  % Frequência de parada abaixo da faixa de passagem (5 MHz)
Fstop2 = 25e6; % Frequência de parada acima da faixa de passagem (25 MHz)
Rp = 1;        % Ripple na faixa de passagem (dB)
Rs = 60;       % Atenuação na faixa de rejeição (dB)

% Projetar o filtro passa-banda
d = designfilt('bandpassiir', ...
               'StopbandFrequency1',Fstop1, ...
               'PassbandFrequency1',Fpass1, ...
               'PassbandFrequency2',Fpass2, ...
               'StopbandFrequency2',Fstop2, ...
               'StopbandAttenuation1',Rs, ...
               'StopbandAttenuation2',Rs, ...
               'PassbandRipple',Rp, ...
               'SampleRate',Fs);

% Visualizar as características do filtro
fvtool(d, 'Fs', Fs); % Resposta em frequência

% Testar o filtro com um sinal de exemplo
t = 0:1/Fs:1e-6; % Tempo de 1 µs
x = cos(2*pi*15e6*t) + cos(2*pi*5e6*t); % Sinal de teste com frequências de 15 MHz e 5 MHz

% Filtrar o sinal
y = filter(d, x);

% Plotar o sinal original e o sinal filtrado
figure;
subplot(2,1,1);
plot(t, x);
title('Sinal Original');
xlabel('Tempo (s)');
ylabel('Amplitude');

subplot(2,1,2);
plot(t, y);
title('Sinal Filtrado');
xlabel('Tempo (s)');
ylabel('Amplitude');
