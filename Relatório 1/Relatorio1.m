%% Gerando Vetores
clear all
close all
clc 
%{
Se discreto = 1 serão os valores discretizados ( comando stem)
Se discreto = 0 será utilizado o comando plot
Se erroAcumuladoPlot = 1 mostra o erroAcumulado em um gráfico
se erroAcumuladoPrint = 1 imprime o valor do erroAcumulado no terminal
%}
discreto = 0;
erroAcumuladoPlot = 1;
erroAcumuladoPrint = 1;

% intervalo de tempo entre amostras
deltaT = 1e-4;

T0 = 0;      % instante inicial
Tf = 0.02;   % instante final
Tempo = T0:deltaT:Tf; % vetor tempo

% vetor de amostras para x(t) = cos 100πt
Xt = cos(100*pi*Tempo);


resol = [2 4 6 8 16]; % número de bits
erroAcumulado = [0 0 0 0 0]; %inicializado vetor de erroAcumulado

figure (1)
subplot(2,1,2)

for i = 1:length(resol)
    mi = min(Xt);
    Xi = Xt - mi; % retirar os valores negativos da curva
    ma = max(Xi);
    Xnorm = Xi ./ ma; % normalizar a curva compatível com o número de níveis
    Yt = ma * round(Xnorm * (2^(resol(i))-1)) / (2^(resol(i))-1); % retira a normalização
    Yt = Yt + round(mi); % retorna os valores negativos
    
    %Faz a soma do erro absoluto para cada resolução diferente.
    erroAcumulado(i) = sum(abs(abs(Xt)-abs(Yt)));
    %Plota os valores de acordo com o método escolhido no inicio
    if (discreto == 1)
    stem(Tempo, Yt)
    else
    plot(Tempo, Yt)
    end
    hold on
end
title("Sinais convertidos");
% Adicionar legenda
hLegend = legend([num2str(resol(1)) ' bits'], [num2str(resol(2)) ' bits'],...
                    [num2str(resol(3)) ' bits'],[num2str(resol(4)) ' bits'],...
                        [num2str(resol(5)) ' bits']);
% Centralizar a legenda na parte superior do gráfico
set(hLegend, 'Position', [0.5, 0.85, 0.2, 0.2], 'Units', 'normalized', 'Location', 'north');
grid on
subplot(2,1,1)
plot(Tempo,Xt);
title("Sinal original");
grid on
if (erroAcumuladoPlot == 1)
figure (2)
stem(resol, erroAcumulado, 'filled');
title("Erro Acumulado para diferentes resoluções")
grid on
% Definir os limites dos eixos
xlim([1.7 17]) % Limite do eixo x de -0.5 até 17
ylim([-0.2 35])   % Limite do eixo y de -1 até 35
end

if(erroAcumuladoPrint == 1)
%Imprime o vetor dos diferentes erros acumulados no terminal
erroAcumulado
end
