%clear all
Fs = 40 * 10^3;  % sample rate 
settings;
% carga al workspace variables de utilidad
archivo_para_transmitir=fopen('input_data.txt');
texto_leido = fread(archivo_para_transmitir);
closedfile=fclose('all');
texto_leido=transpose(texto_leido);
texto_leido=double(texto_leido);
texto_leido=dec2bin(texto_leido);
[row, col , deep]= size(texto_leido);  % dimension imagen
signal = []; % mensaje total a enviar
%% ********* HEADER *********
	% Secuencia de 5k, 6k, y 7k 
header_freq=[header_freq1 header_freq2 header_freq3];
head_t = 0:1/Fs:head_dur; 
header = [sin(2*pi*header_freq(1)*head_t), sin(2*pi*header_freq(2)*head_t), sin(2*pi*header_freq(3)*head_t)];
signal = [signal, header];
% tiempo2=0:1/Fs:(length(header)-1)/Fs;
% figure
% plot(tiempo2,header);

%% ********* Codificacion de dimension de la imagen *********
tam_t = 0:1/Fs:tam_dur;
header = sin(2*pi*(5*row+4650)*tam_t);
signal = [signal, header];
% figure
% plot(tam_t,header);
%% ********* MESSAGE *********
% Acomodar datos para emitir
vector_binario_datos=zeros(1,row*col);
for counter=1:(row)
    for counter2=1:col
        vector_binario_datos(1,(counter-1)*7+counter2)=str2double(texto_leido(counter,counter2));
    end
end

%% "Modular"
t=0:1/Fs:pix_dur;
for counter=1:length(vector_binario_datos)
    if vector_binario_datos(counter)==1
        frecuencia=freq1;
    else
        frecuencia=freq2;
    end
    sample(1,(counter-1)*length(t)+1:counter*length(t))=sin(frecuencia*2*pi.*t);
end
signal = [signal, sample];
% Canal ideal
disp('presiona a una tecla para continuar')
% save('audio.mat','signal');
 pause
 sound(signal, Fs)


L = length(signal);
NFFT = 2^nextpow2(L);
Y = fft(signal, NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
plot(f, 2*abs(Y(1:NFFT/2+1)));
xlabel('Frecuencia [Hz]')
ylabel('Amplitud')