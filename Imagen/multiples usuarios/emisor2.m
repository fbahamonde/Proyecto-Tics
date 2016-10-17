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
        pix=[img(i,j,2)];
        pix = double(pix);
		frecP = [p2f(pix,GbaseF,dF)];
		sample = sin(2 * pi * frecP* pix_t);
		%sample = sum(sample);
  		signal = [signal, sample];
    end
end
%save('audiog.mat','signal');
disp('presiona a una tecla para continuar')
pause
sound(signal, Fs)
%wavwrite(signal,Fs,'audio.wav')

