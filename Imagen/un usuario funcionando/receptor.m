close all
clear
clc
%{
Receptor    
recibe el sonido y decodifica la imagen
-graba por un tiempo
-detecta el header
-decodifica dimension de la imagen
-recupera los pixeles
%}
% ********* Grabacion *********
fs = 40 * 10^3;  % sample rate
disp('Recibiendo Se�al')
t_f = 95; % duracion de la grabacion (segs)
recobj = audiorecorder(fs, 16, 1); 
recordblocking(recobj, t_f);
signal = recobj.getaudiodata;
% load('audio.mat');  % prueba canal perfecto
settings;
% carga al workspace:
    % head_dur
    % tam_dur
    % pix_dur
    % RbaseF
    % GbaseF
    % BbaseF
    % dF
% ********* Header *********
head_t = 0:1/fs:head_dur; 
true_head = [sin(2*pi*5000*head_t), sin(2*pi*6000*head_t), sin(2*pi*7000*head_t)];
test_dur = 4; % segs
test_head = signal(1:test_dur*fs); % header a detectar (4 segundos de ventana)
% 'acor' contiene la correlacion para cada desplazamiento de 'lag'
[acor, lag] = xcorr(test_head, true_head); % correlacion cruzada
[~,I] = max(abs(acor));  % maxima correlacion
lagDiff = lag(I);  % desviacion en frames donde esta la maxima corr

% ********* Decodificacion dimension imagen *********
% head_init tiene el inicio de donde empieza a codificarse la dimension
head_init = abs(lagDiff) + length(true_head);  
L = tam_dur * fs;  % frames que dura la info de la dimension
body = signal(head_init:head_init + L -1); % info dimension
Z = fft(body);
Z1 = Z(1:(L/2)+1);
Z1(2:end-1) = 2*Z1(2:end-1);
Z1 = abs(Z1);
Z12 = abs(Z1);
efe = fs*(0:(L/2))/L;
index = (10000 > efe) & (efe >= 5600); % indices filtrado
efe = efe(index);
Z1 = Z1(index); % filtrado ;)
[~,f_i] = max(Z1);
efe(f_i);  % encuentra el tono

efe2 = fs*(0:(L/2))/L;
index2 = (3000> efe2) & (efe2 > 2000); % indices filtrado
efe2 = efe2(index2);
Z12 = Z12(index2); % filtrado ;)
[~,f_i2] = max(Z12);
efe2(f_i2)  % encuentra el tono

tam = round((efe(f_i)-5600)/5)
tam2 = round((efe2(f_i2)-2000)/5)
%hasta aqui funcionando
% ********* Mensaje *********
init_time = head_init + length(body);
red = zeros(1, tam*tam2);
green = zeros(1, tam*tam2);
blue = zeros(1, tam*tam2);
pix_t=0:1/fs:pix_dur; 
pix_wid = length(pix_t);
f_val = (0:pix_wid/2) * fs/pix_wid;  % frecuencias de la fft
indexR =((RbaseF+255*dF+10)> f_val) & (f_val > RbaseF);
indexG =((GbaseF+255*dF+10)> f_val) & (f_val > GbaseF);
indexB =((BbaseF+255*dF+10)> f_val) & (f_val > BbaseF);

Rf_val = f_val(indexR);
Gf_val = f_val(indexG);
Bf_val = f_val(indexB);

message = signal(init_time:end);
for i = 1:(tam*tam2)
        sample = message((i-1)*pix_wid+1:i*pix_wid);
        win = hamming(pix_wid);
        sample = sample.*win;
        Y = fft(sample); % saca fft por columnas
        Y1 = Y(1:round(pix_wid/2) + 1);        
        Y1 = abs(Y1);
        Y1r = Y1(indexR);
        Y1r = abs(Y1r); % corta banda del rojo
        [~, indr] = max (Y1r); % encuentra máxima frecuencia en la banda del rojo
        red(i) = f2p(Rf_val(indr),RbaseF, dF); % guarda el valor r del pixel
        Y1g= Y1(indexG);
        Y1g = abs(Y1g); % corta banda del verde
        [frecg,indg] = max (Y1g); % encuentra máx frec en la banda del verde
        green(i) = f2p(Gf_val(indg),GbaseF, dF); % guarda el valor g del pixel
        Y1b= Y1(indexB); 
        Y1b = abs(Y1b); % corta banda del azul
        [frecb,indb] = max (Y1b); % encuentra máx frec de la banda del azul
        blue(i) = f2p(Bf_val(indb),BbaseF, dF); % guarda el valor b del pixel
end
red8 = uint8(reshape(red, tam, tam2)); % Matriz con valores de rojo
green8 = uint8(reshape(green, tam, tam2)); % Matriz con valores de verde
blue8 = uint8(reshape(blue, tam, tam2)); % Matriz con valores de azul

imgrec = cat(3,red8',green8',blue8'); % Imagen reconstruida
% true_img = imread(foto);
% subplot(2,1,1)
% imshow(imgrec);
% xlabel('Recibida')
% subplot(2,1,2)
% imshow(true_img);
% xlabel('Original')
imshow(imgrec)


