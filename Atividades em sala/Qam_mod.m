function y = Qam_mod(x,fc,fs,opc)
ts = 1/fs;
Ns = length(x);
t=0:ts:ts*(Ns-1);
y = x.*cos(2*pi*fc.*t) + opc.*sin(2*pi*fc.*t);
end