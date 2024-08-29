clc
clear
A=1;
nsb=16;
Nb=10000;
sb=randi([0,1],1,Nb); %sequência de bits
gt=A*sb;            %sinal
rt= randn(1,10000); %ruído gaussiano
xt = gt+rt;         %sinal recebido

ps = mean(gt.^2); %potência do sinal = A^2/2
pn = mean(rt.^2); %potência do ruído

SNR = 10*log10(ps/pn); %Potência do sinal em Db

%Calculando A, considerando que a potência de ruído não muda
ps2 = zeros(1,12);
for SNR2 = -3:1:10;  %Dbs desejados
    ps2(SNR2+4)= (10^(SNR2/10))*pn;  %Potências do sinal necessárias
end
A2 = reshape(((ps2*2).^(1/2)),14,1); %Amplitudes a partir das potências

%Sinais com amplitudes que variam a SNR de -3 a 10 dB comparadas ao ruído
lambda=0;
gt1=A2(1).*sb;
xt1 = (gt1+rt)>lambda;         %sinal recebido (em binário)

gt2=A2(2).*sb;
xt2 = (gt2+rt)>lambda;         %sinal recebido

gt3=A2(3).*sb;
xt3 = (gt3+rt)>lambda;         %sinal recebido

gt4=A2(4).*sb;
xt4 = (gt4+rt)>lambda;         %sinal recebido

gt5=A2(5).*sb;
xt5 = (gt5+rt)>lambda;         %sinal recebido

gt6=A2(6).*sb;
xt6 = (gt6+rt)>lambda;         %sinal recebido

gt7=A2(7).*sb;
xt7 = (gt7+rt)>lambda;         %sinal recebido

gt8=A2(8).*sb;
xt8 = (gt8+rt)>lambda;         %sinal recebido

gt9=A2(9).*sb;
xt9 = (gt9+rt)>lambda;         %sinal recebido

gt10=A2(10).*sb;
xt10 = (gt10+rt)>lambda;         %sinal recebido

gt11=A2(11).*sb;
xt11 = (gt11+rt)>lambda;         %sinal recebido

gt12=A2(12).*sb;
xt12 = (gt12+rt)>lambda;         %sinal recebido

gt13=A2(13).*sb;
xt13 = (gt13+rt)>lambda;         %sinal recebido

gt14=A2(14).*sb;
xt14 = (gt14+rt)>lambda;         %sinal recebido

BERm3= biterr(xt1,sb); %Biterror -3dB
BERm2= biterr(xt2,sb);
BERm1= biterr(xt3,sb);
BER0= biterr(xt4,sb); %Biterror 0dB
BER1= biterr(xt5,sb);
BER2= biterr(xt6,sb);
BER3= biterr(xt7,sb);
BER4= biterr(xt8,sb);
BER5= biterr(xt9,sb);
BER6= biterr(xt10,sb);
BER7= biterr(xt11,sb);
BER8= biterr(xt12,sb);
BER9= biterr(xt13,sb);
BER10= biterr(xt14,sb);

x=[BERm3, BERm2, BERm1, BER0, BER1, BER2, BER3, BER4, BER5, BER6, BER7, BER8, BER9, BER10];
plot(-3:1:10,x) %SNR por Db



%Parte 2
bit1 = repmat([ones(1,nsb/2) -1*ones(1,nsb/2)],1,Nb); %Bit1
bit0 = repmat([-1*ones(1,nsb/4) ones(1, nsb/2) -1*ones(1,nsb/4)],1,Nb); %Bit0
fc1 = [-ones(1,nsb/2) 1*ones(1,nsb/2)]; %Filtro casado do 1
fc0 = [-1*ones(1,nsb/4) ones(1, nsb/2) -1*ones(1,nsb/4)];   %Filtro casado do 0
snrzu1 = (reshape(repmat(sb, nsb,1),1,nsb*Nb)); %sinal nao retorno 0, bit 1
snrzu0 = -snrzu1+1; %sinal nao retorno 0, bit 0



ruidort = (reshape(repmat(rt, nsb,1),1,nsb*Nb)); %ruido 1x16000 amostras
sinalgt = A*(snrzu1.*bit1 + snrzu0.*bit0); %sinal 1x16000 amostras

sinalrxt = ruidort+sinalgt; %sinal recebido

ps2 = mean(sinalgt.^2); %potência do sinal
pn2 = mean(ruidort.^2); %potência do ruído
SNR2 = 10*log10(ps2/pn2);


%Calculando A, considerando que a potência de ruído não muda
for SNR3 = -3:1:10;  %Dbs desejados
    ps3(SNR3+4)= (10^(SNR3/10))*pn2;  %Potências do sinal necessárias
end
A3 = reshape(((ps3).^(1/2)),14,1); %Amplitudes a partir das potências



sinalgt1 = A3(1)*(snrzu1.*bit1 + snrzu0.*bit0); %sinal 1x16000 amostras
sinalrxt1 = sinalgt1+ruidort;

sinalgt2 = A3(2)*(snrzu1.*bit1 + snrzu0.*bit0); %sinal 1x16000 amostras
sinalrxt2 = sinalgt2+ruidort;

sinalgt3 = A3(3)*(snrzu1.*bit1 + snrzu0.*bit0); %sinal 1x16000 amostras
sinalrxt3 = sinalgt3+ruidort;

sinalgt4 = A3(4)*(snrzu1.*bit1 + snrzu0.*bit0); %sinal 1x16000 amostras
sinalrxt4 = sinalgt4+ruidort;

sinalgt5 = A3(5)*(snrzu1.*bit1 + snrzu0.*bit0); %sinal 1x16000 amostras
sinalrxt5 = sinalgt5+ruidort;

sinalgt6 = A3(6)*(snrzu1.*bit1 + snrzu0.*bit0); %sinal 1x16000 amostras
sinalrxt6 = sinalgt6+ruidort;

sinalgt7 = A3(7)*(snrzu1.*bit1 + snrzu0.*bit0); %sinal 1x16000 amostras
sinalrxt7 = sinalgt7+ruidort;

sinalgt8 = A3(8)*(snrzu1.*bit1 + snrzu0.*bit0); %sinal 1x16000 amostras
sinalrxt8 = sinalgt8+ruidort;

sinalgt9 = A3(9)*(snrzu1.*bit1 + snrzu0.*bit0); %sinal 1x16000 amostras
sinalrxt9 = sinalgt9+ruidort;

sinalgt10 = A3(10)*(snrzu1.*bit1 + snrzu0.*bit0); %sinal 1x16000 amostras
sinalrxt10 = sinalgt10+ruidort;

sinalgt11 = A3(11)*(snrzu1.*bit1 + snrzu0.*bit0); %sinal 1x16000 amostras
sinalrxt11 = sinalgt11+ruidort;

sinalgt12 = A3(12)*(snrzu1.*bit1 + snrzu0.*bit0); %sinal 1x16000 amostras
sinalrxt12 = sinalgt12+ruidort;

sinalgt13 = A3(13)*(snrzu1.*bit1 + snrzu0.*bit0); %sinal 1x16000 amostras
sinalrxt13 = sinalgt13+ruidort;

sinalgt14 = A3(14)*(snrzu1.*bit1 + snrzu0.*bit0); %sinal 1x16000 amostras
sinalrxt14 = sinalgt14+ruidort;


srconv = conv(sinalrxt1,fc1); %A cada T amostras (16 já que nsb=16), ou tem um valor muito próximo de 0, ou muito próximo de 1
srconv0 = conv(sinalrxt14,fc0);

filtrado=srconv(nsb:nsb:nsb*Nb)>0.5;

erro=biterr(filtrado,sb);