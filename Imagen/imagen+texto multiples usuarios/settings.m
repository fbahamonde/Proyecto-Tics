% duracion de la grabacion (segs)
t_f = 100; 
% duracion (tiempo) de cada tono del header
head_dur = 0.6; 
% duracion (tiempo) donde esta el dimension de la imagen
tam_dur = 0.3;
% duracion (tiempo) de ventana donde esta el pixel
pix_dur = 0.35; % segs 
% funciones para pasar de frecuencia al valor del pixel 
RbaseF = 3000;
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
ftr1=14000;
ftr2=14500;
ftg1=12000;
ftg2=12500;
ftb1=13000;
ftb2=13500;
foto='angel.png';

%% Frecuencias para emitir texto
freq1=2000;
freq2=2500;
f1r=5000/2;%5000
f2r=3500; %6k
f3r=4500; %7k
header_freq1=f1r;
header_freq2=f2r;
header_freq3=f3r;
%frecuencias para el tamaño
ftr1=4650; %14k
