sidqamrx = sinalqamrx.*cos(2*pi*fc.*t) - 1i*sinalqamrx.*sin(2*pi*fc.*t);
sidpskrx = sinalpskrx.*cos(2*pi*fc.*t) - 1i*sinalpskrx.*sin(2*pi*fc.*t);
spskrx = resample(sidpskrx,Q,P,10);
sqamrx = resample(sidqamrx,Q,P,10);