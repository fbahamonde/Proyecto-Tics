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
%frecuencias para los header
f1r=5000;
f2r=6000;
f3r=7000;
f1g=2000;
f2g=3000;
f3g=4000;
f1b=8000;
f2b=9000;
f3b=10000;
%frecuencias para el tamaño
ftr1=5600;
ftr2=2000;
ftg1=12000;
ftg2=12500;
ftb1=13000;
ftb2=13500;
foto='cara.jpg';