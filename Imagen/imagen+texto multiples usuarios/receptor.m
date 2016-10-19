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
                        %% settings
                        settings;
                        %{
% carga al workspace:
                        % head_dur
                        % tam_dur
                        % pix_dur
                        % RbaseF
                        % GbaseF
                        % BbaseF
                        % dF
%}
%%                    ********* Grabacion *********
                        fs = 40 * 10^3;  % sample rate
                        disp('Recibiendo señal')
                        disp(['Duración Grabación :' num2str(t_f) ' seg.'])                    

                        recobj = audiorecorder(fs, 16, 1);
                        recordblocking(recobj, t_f);
                        signal = recobj.getaudiodata;
%                         %% Canal ideal
%                         load('audiob3.mat');
%                         load('audiog.mat');
%                         load('audio.mat');% prueba canal perfecto
%                         if length(signal)<length(signal2)
%                             parche=zeros(1,length(signal2)-length(signal));
%                             signal=[signal,parche];
%                             length(signal);
%                         end
%                         signal=signal2+signal3+signal;  
                        

%%                                  Receptor texto
% Header
head_t = 0:1/fs:head_dur; 
true_head = [sin(2*pi*header_freq1*head_t), sin(2*pi*header_freq2*head_t), sin(2*pi*header_freq3*head_t)];
test_dur = 4; % segs
test_head = signal(1:test_dur*fs); % header a detectar (4 segundos de ventana)
% 'acor' contiene la correlacion para cada desplazamiento de 'lag'
[acor, lag] = xcorr(true_head, test_head); % correlacion cruzada
[~,I] = max(abs(acor));  % maxima correlacion
lagDiff = lag(I);  % desviacion en frames donde esta la maxima corr
% Decodificacion dimension imagen
% head_init tiene el inicio de donde empieza a codificarse la dimension
head_init = abs(lagDiff) + length(true_head);  
L = round(tam_dur * fs);  % frames que dura la info de la dimension
body = signal(head_init:head_init + L -1); % info dimension
Z = fft(body);
Z1 = Z(1:round(L/2) -1);
Z1 = abs(Z1);
efe =  (0:round(L/2)) * fs/L;
index = (4800 > efe) & (efe > 4650); % indices filtrado
efe = efe(index);
Z1 = Z1(index); % filtrado ;)
[~,f_i] = max(Z1);
efe(f_i); % encuentra el tono
tam = round((efe(f_i) - ftr1)/5 )*7; %Largo del vector con informacion
disp(['Cantidad de caracteres del texto: ' num2str(tam/7)]);


%  Mensaje
init_time = head_init + length(body);
datos_recibidos=signal((init_time+1):(init_time+(pix_dur*fs*tam)));
largo_ventana=fs*pix_dur;
reconstruccion_mensaje=zeros(1,tam);
for contador2=1:tam
    datos=datos_recibidos(((contador2-1)*largo_ventana)+1:(largo_ventana*contador2));
    L = length(datos);
    NFFT = 2^nextpow2(L);
    Y = fft(datos, NFFT)/L;
    f = fs/2*linspace(0,1,NFFT/2+1);
    index = (3100 > f) & (f > 1900); % indices filtrado
    Y=Y(index);
    f=f(index);
    [~, indicefreq]=max(abs(Y));
    frecuencia_maxima_ventana=f(indicefreq);
    if abs(frecuencia_maxima_ventana-freq1)<50
        reconstruccion_mensaje(1,contador2)=1;
    else
        reconstruccion_mensaje(1,contador2)=0;
    end
end
probando=datos_recibidos(1:end);
L = length(probando);
NFFT = 2^nextpow2(L);
Y = fft(probando, NFFT)/L;
f = fs/2*linspace(0,1,NFFT/2+1);
figure
plot(f, 2*abs(Y(1:NFFT/2+1)));
%%agrupar de a 7
counter_2=1;
salida=zeros(tam/7,7);
w=reconstruccion_mensaje;
while 1
    if (counter_2*7)>length(w)
        break
    end
    salida(counter_2,:)=w((counter_2-1)*7+1:(counter_2*7));
    counter_2=counter_2+1;
end

%%
                                    %Matriz Green and Red
% ********* Header *********
head_t = 0:1/fs:head_dur; 
true_head = [sin(2*pi*f1g*head_t), sin(2*pi*f2g*head_t), sin(2*pi*f3g*head_t)];
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
index = (ftg1+400 > efe) & (efe > ftg1); % indices filtrado
efe = efe(index);
Z1 = Z1(index); % filtrado ;)
[~,f_i] = max(Z1);
efe(f_i);  % encuentra el tono
efe2 = fs*(0:(L/2))/L;
index2 = (ftg2+400> efe2) & (efe2 > ftg2); % indices filtrado
efe2 = efe2(index2);
Z12 = Z12(index2); % filtrado ;)
[~,f_i2] = max(Z12);
efe2(f_i2);  % encuentra el tono
colr = round((efe(f_i)-ftg1)/5);
filr = round((efe2(f_i2)-ftg2)/5);
disp(['Dimensión matriz R: ' num2str(colr) 'x' num2str(filr)]);
disp(['Dimensión matriz G: ' num2str(colr) 'x' num2str(filr)]);
% ********* Mensaje *********
init_time = head_init + length(body);
green = zeros(1, colr*filr);
pix_t=0:1/fs:pix_dur; 
pix_wid = length(pix_t);
f_val = (0:pix_wid/2) * fs/pix_wid;  % frecuencias de la fft
indexG =((GbaseF+255*dF+10)> f_val) & (f_val > GbaseF);

Gf_val = f_val(indexG);

message = signal(init_time:end);
for i = 1:(colr*filr)
        sample = message((i-1)*pix_wid+1:i*pix_wid);
        win = hamming(pix_wid);
        sample = sample.*win;
        Y = fft(sample); % saca fft por columnas
        Y1 = Y(1:round(pix_wid/2) + 1);        
        Y1 = abs(Y1);
        Y1g= Y1(indexG);
        Y1g = abs(Y1g); % corta banda del verde
        [frecg,indg] = max (Y1g); % encuentra mÃ¡x frec en la banda del verde
        green(i) = f2p(Gf_val(indg),GbaseF, dF); % guarda el valor g del pixel
end
green8 = uint8(reshape(green, colr, filr)); % Matriz con valores de verde

% ********* Mensaje *********
init_time = head_init + length(body);
red = zeros(1, colr*filr);
pix_t=0:1/fs:pix_dur; 
pix_wid = length(pix_t);
f_val = (0:pix_wid/2) * fs/pix_wid;  % frecuencias de la fft
indexR =((RbaseF+255*dF+10)> f_val) & (f_val > RbaseF);
Rf_val = f_val(indexR);
message = signal(init_time:end);
for i = 1:(colr*filr)
        sample = message((i-1)*pix_wid+1:i*pix_wid);
        win = hamming(pix_wid);
        sample = sample.*win;
        Y = fft(sample); % saca fft por columnas
        Y1 = Y(1:round(pix_wid/2) + 1);        
        Y1 = abs(Y1);
        Y1r = Y1(indexR);
        Y1r = abs(Y1r); % corta banda del rojo
        [~, indr] = max (Y1r); % encuentra mÃ¡xima frecuencia en la banda del rojo
        red(i) = f2p(Rf_val(indr),RbaseF, dF); % guarda el valor r del pixel
end
red8 = uint8(reshape(red, colr, filr)); % Matriz con valores de rojo
%%
                                        %Matriz Blue
% ********* Header *********
head_t = 0:1/fs:head_dur; 
true_head = [sin(2*pi*f1b*head_t), sin(2*pi*f2b*head_t), sin(2*pi*f3b*head_t)];
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
index = (ftb1+400 > efe) & (efe > ftb1); % indices filtrado
efe = efe(index);
Z1 = Z1(index); % filtrado ;)
[~,f_i] = max(Z1);
efe(f_i);  % encuentra el tono
efe2 = fs*(0:(L/2))/L;
index2 = (ftb2+400> efe2) & (efe2 > ftb2); % indices filtrado
efe2 = efe2(index2);
Z12 = Z12(index2); % filtrado ;)
[~,f_i2] = max(Z12);
efe2(f_i2);  % encuentra el tono
colb = round((efe(f_i)-ftb1)/5);
filb = round((efe2(f_i2)-ftb2)/5);
disp(['Dimensión matriz B: ' num2str(colb) 'x' num2str(filb)]);

% ********* Mensaje *********
init_time = head_init + length(body);
blue = zeros(1, colb*filb);
pix_t=0:1/fs:pix_dur; 
pix_wid = length(pix_t);
f_val = (0:pix_wid/2) * fs/pix_wid;  % frecuencias de la fft
indexB =((BbaseF+255*dF+10)> f_val) & (f_val > BbaseF);
Bf_val = f_val(indexB);
message = signal(init_time:end);
for i = 1:(colb*filb)
        sample = message((i-1)*pix_wid+1:i*pix_wid);
        win = hamming(pix_wid);
        sample = sample.*win;
        Y = fft(sample); % saca fft por columnas
        Y1 = Y(1:round(pix_wid/2) + 1);        
        Y1 = abs(Y1);
        Y1b= Y1(indexB); 
        Y1b = abs(Y1b); % corta banda del azul
        [frecb,indb] = max (Y1b); % encuentra mÃ¡x frec de la banda del azul
        blue(i) = f2p(Bf_val(indb),BbaseF, dF); % guarda el valor b del pixel
end
blue8 = uint8(reshape(blue, colb, filb)); % Matriz con valores de azul
%% Resultados
disp(['dimension de la imagen: ' num2str(filb) 'x' num2str(colb)])
imgrec = cat(3,red8',green8',blue8'); % Imagen reconstruida
true_img = imread(foto);
figure (1)
subplot(1,2,1)
imshow(imgrec);
xlabel('Recibida')
subplot(1,2,2)
imshow(true_img);
xlabel('Original')

mensaje=transpose(char(bin2dec(num2str(salida))))
disp(['Mensaje: ' mensaje]);

fileID = fopen('texto_recibido.txt','w');
fprintf(fileID, mensaje);
fclose(fileID);

L = length(signal);
NFFT = 2^nextpow2(L);
Y = fft(signal, NFFT)/L;
f = fs/2*linspace(0,1,NFFT/2+1);
figure (2)
plot(f, 2*abs(Y(1:NFFT/2+1)));
xlabel('Frecuencia [Hz]')
ylabel('Amplitud')

fileID = fopen('texto_recibido.txt','w');
fprintf(fileID, mensaje);
fclose(fileID);

