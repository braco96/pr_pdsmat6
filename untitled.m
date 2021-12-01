[x,fs]=audioread("digitos.wav");

eje_t = 1/fs:1/fs:length(x)/fs;

figure(1);
plot(eje_t,x);

x_fil = x - mean(x);

figure(1);
hold on;
plot(eje_t, x_fil);
hold off;

n_w = 512;
n_p = 128;
solapa = 256;
desplaza = n_w - solapa;


figure(2);
spectrogram(x, n_w, solapa, n_w, 'yaxis', fs);
colormap 'jet';
colorbar
figure(3);
spectrogram(x_fil, n_w, solapa, n_w, 'yaxis', fs);
colormap 'jet';

ini = 14500;
final = ini + n_w -1;
trozo = x_fil(ini:final);
final1 = ini + n_p -1;
trozop = x_fil(ini:final1);

figure(4);
plot(trozo);
trozo_w = trozo.*hamming(n_w);

figure(4);
hold on;
plot(trozo_w, 'r');
hold off;

% Mayus por ser frecuencia
ceros = zeros(1,n_w)
TROZO = fft(trozo);
TROZO_DB = 20*log10(abs(TROZO));
TROZO_W = fft(trozo_w);
TROZO_W_DB = 20*log10(abs(TROZO_W));

figure(5);
plot(TROZO_DB);
hold on;
plot(TROZO_W_DB,'r');
hold off;

% Quedaba por calcular trozo_p, no he podido copiarlo
trozo_p = [trozo',ceros,ceros,ceros]';
trozo_w_p = [trozo_w',ceros,ceros,ceros]';

TROZO_P = fft(trozo_p);
TROZO_P_DB = 20*log10(abs(TROZO_P));
TROZO_W_P = fft(trozo_w_p);
TROZO_W_P_DB = 20*log10(abs(TROZO_W_P));
figure(6);

hold on;
plot(TROZO_P_DB,'r');
hold off;

descarta=n_w/desplaza-1;
n_ventanas=floor(length(x)/desplaza)-descarta;
matriz_espec=zeros(n_w/2+1,n_ventanas);
n_cruces=zeros(1,n_ventanas);
for i=1 : n_ventanas
    %Cortamos la señal en trozos.
    ini=(i-1)*desplaza+1;
    final=ini+n_w-1;
    trozo=x(ini:final);
    trozo_f=x_fil(ini:final);
    %2 Enventanado.
    trozo_f_w=trozo_f.*hamming(n_w);
     %3 Procesado.
      %3.1 Calculo del espectograma.
      TROZO_F_W=abs(fft(trozo_f_w));
      TROZO_F_W_DB=20*log10(TROZO_F_W)
      matriz_espec(:,i)=TROZO_F_W_DB(1:n_w/2+1);
      %3.2Calculo de cruces por cero.
      for j =2:n_w
          if(trozo_f(j)>=0 & trozo_f(j-1)<0)
              n_cruces(i)=n_cruces(i)+1;
          end
      end
      %3.3 Calculo de la energia :
      energias(i)=sum(trozo_f.*trozo_f);
end

figure(7);
surf(matriz_espec);
shading interp;
colormap 'jet';
figure(8);
plot(n_cruces)
figure(9);
plot(energias);


