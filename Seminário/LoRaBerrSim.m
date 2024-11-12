close all
clear all
clc

SF = 7:1:12;
SNR = -20;
BW = 1000; % largura de banda
Fs = 1000; %amostragem
s = 100; %send symbol '100'
%--- Gera o símbolo
for b=1:6
    num_samples = (2^SF(b))*Fs/BW;
    k = s; %add s to k to start (defines the data symbol)
    lora_symbol = zeros(1,num_samples);
    for n=1:num_samples
        if k>=(2^SF(b)) % Este if cumpre a função do operador mod
            k = k-2^SF(b);
        end
        k=k+1;
        lora_symbol(n) = (1/(sqrt(2^SF(b))))*exp(1i*2*pi*(k)*(k/(2^SF(b)*2)));
    end
    for j=1:100
        %add noise
        lora_symbol_noisy = awgn(lora_symbol, SNR, 'measured');
        %transmit
        %Receiver below
        %---- Generate the Base Down Chirp
        base_down_chirp = zeros(1,num_samples);
        k=0;
        for n=1:num_samples
            if k>= (2^SF(b))
                k = k-2^SF(b);
            end
            k = k+1;
            base_down_chirp(n) = (1/(sqrt(2^SF(b))))*exp(-1i*2*pi*(k)*(k/(2^SF(b)*2)));
        end
        dechirped = lora_symbol_noisy .* (base_down_chirp);
        corrs = (abs(fft(dechirped)).^2);
%         plot(corrs)
        [~,ind] = max(corrs);
        ind2(j) = ind;
%         pause(0.01)
    end
symbol_error_rate(b) = sum(ind2~=s+1)/j    
end

stem(SF, symbol_error_rate)
xlim([5.8 12.3])
grid on

