% Especificações do filtro
Fs = 100e6; % Frequência de amostragem (100 MHz)
Fpass1 = 10e6; % Frequência de passagem baixa (10 MHz)
Fpass2 = 20e6; % Frequência de passagem alta (20 MHz)

% Frequências angulares
omega1 = 2 * pi * Fpass1;
omega2 = 2 * pi * Fpass2;

% Coefficients for a bandpass filter in analog domain
% Standard bandpass filter coefficients
% s^2 + omega1^2 in numerator and s^2 + omega2^2 in denominator

% Define the poles and zeros
num = [1 0 -omega1^2]; % Numerator polynomial coefficients (s^2 + omega1^2)
den = [1 0 -omega2^2]; % Denominator polynomial coefficients (s^2 + omega2^2)

% Display poles and zeros
[z, p, k] = tf2zp(num, den);

disp('Zeros:');
disp(z);
disp('Polos:');
disp(p);

% Plot the poles and zeros
figure;
zplane(z, p);
title('Polos e Zeros do Filtro Passa-Banda');
xlabel('Parte Real');
ylabel('Parte Imaginária');
