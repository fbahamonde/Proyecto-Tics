close all
clear
clc
%{
Lee imagen
emite linea por linea pixel por pixel de izquierda a derecha de arriba a abajo
%}
%********* PARAMETROS *********
Fs = 40 * 10^3;  % sample rate 
settings;
% carga al workspace:
	% head_dur
	% tam_dur
	% pix_dur
	% RbaseF
	% GbaseF
	% BbaseF
	% dF
img = imread(foto);
[row, col , deep]= size(img);  % dimension imagen
signal = []; % mensaje total a enviar

%********* HEADER *********
	% Secuencia de 5k, 6k, y 7k 
head_t = 0:1/Fs:head_dur; 
header = [sin(2*pi*f1g*head_t), sin(2*pi*f2g*head_t), sin(2*pi*f3g*head_t)];
signal = [signal, header];

%********* Codificacion de dimension de la imagen *********
tam_t = 0:1/Fs:tam_dur;
header = sin(2*pi*(ftg1+col*5)*tam_t)+ sin(2*pi*(ftg2+row*5)*tam_t);
g=length(header);
ff=fft(header);
Z1 = ff(1:(g/2)+1);
Z1(2:end-1) = 2*Z1(2:end-1);
Z1 = abs(Z1);
signal = [signal, header];

%********* MESSAGE *********
pix_t = 0:1/Fs:pix_dur;
for i=1:row
   	for j=1:col
        pix=[img(i,j,1) img(i,j,2)];
        pix = double(pix);
		frecP = [p2f(pix(1),RbaseF,dF), p2f(pix(2),GbaseF,dF)];
		sample = sin(2 * pi * frecP.' * pix_t);
		sample = sum(sample);
  		signal = [signal, sample];
    end
end
signal2=signal;
% save('audiog.mat','signal2');
disp('presiona a una tecla para continuar')
pause
sound(signal, Fs)


signal = [signal, sample];
L = length(signal);
NFFT = 2^nextpow2(L);
Y = fft(signal, NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
plot(f, 2*abs(Y(1:NFFT/2+1)));
xlabel('Frecuencia [Hz]')
ylabel('Amplitud')

