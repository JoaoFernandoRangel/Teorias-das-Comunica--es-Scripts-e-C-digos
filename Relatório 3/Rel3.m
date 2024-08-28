%comentário teste
clear
A=1;
Nb=10000;
sb=randi([0,1],1,Nb); %sequencia de bits
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
A2 = reshape(((ps2*2).^(1/2)),14,1); %Amplitudes à partir das potências

%Sinais com amplitudes que variam a SNR de -3 a 10 dB
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
plot(-3:1:10,x)

