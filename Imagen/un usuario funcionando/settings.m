% duracion (tiempo) de cada tono del header
head_dur = 0.6;
% duracion (tiempo) donde esta el dimension de la imagen
tam_dur = 0.3;
% duracion (tiempo) de ventana donde esta el pixel
pix_dur = 0.35; % segs 
% funciones para pasar de frecuencia al valor del pixel 
RbaseF = 2000;
GbaseF = 6000;
BbaseF = 10000;
dF = 10; % diferencia entre cada valor del pixel
p2f = @(pix, base, df) base + pix*df;
f2p = @(frec, bf, df) round((frec-bf)/df);
foto='angel.png'

