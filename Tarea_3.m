%{
Autor: Antonio Maldonado Pinzon
Tarea 3: Estimacion de Movimiento
Procesamiento de video Ago-Dic 2019
Fecha: 18/09/2019
Programa que calcula el desplazamiento de un objeto en una imagen
%}
clear all;close all;clc;

I1=double(rgb2gray(imread('1.png')));
I21=double(rgb2gray(imread('21.png')));
I22=double(rgb2gray(imread('22.png')));
I23=double(rgb2gray(imread('23.png')));
I24=double(rgb2gray(imread('24.png')));

%******************* FUERZA BRUTA *********************
fprintf('Fuerza Bruta\n');
[dx1 dy1]=FuerzaBruta(I1,I21,11); %Calculamos el desplazamiento de I2_1
fprintf('\nImagen 2_1\nDesplazamiento en filas: %i',dy1)
fprintf('\nDesplazamiento en columnas: %i',dx1)

[dx2 dy2]=FuerzaBruta(I1,I22,11); %Calculamos el desplazamiento de I2_2
fprintf('\n\nImagen 2_2\nDesplazamiento en filas: %i',dy2)
fprintf('\nDesplazamiento en columnas: %i',dx2)

[dx3 dy3]=FuerzaBruta(I1,I23,31); %Calculamos el desplazamiento de I2_3
fprintf('\n\nImagen 2_3\nDesplazamiento en filas: %i',dy3)
fprintf('\nDesplazamiento en columnas: %i',dx3)

[dx4 dy4]=FuerzaBruta(I1,I24,19); %Calculamos el desplazamiento de I2_4
fprintf('\n\nImagen 2_4\nDesplazamiento en filas: %i',dy4)
fprintf('\nDesplazamiento en columnas: %i\n',dx4)
%I1_new=Interpolar(I1,10.3,15.3);

%******************* LUCAS-KANADE **************************
k=1; %Si k=1 se mostrara la animacion del algoritmo
fprintf('\n\nLucas-Kanade\n')
[dx11 dy11]=LucasKanade(I1,I21,k,1); %Calculamos el desplazamiento de I2_1
fprintf('\nImagen 2_1\nDesplazamiento en filas: %f',dy11)
fprintf('\nDesplazamiento en columnas: %f',dx11)

[dx12 dy12]=LucasKanade(I1,I22,k,2); %Calculamos el desplazamiento de I2_2
fprintf('\n\nImagen 2_2\nDesplazamiento en filas: %f',dy12)
fprintf('\nDesplazamiento en columnas: %f',dx12)

[dx13 dy13]=LucasKanade(I1,I23,k,3); %Calculamos el desplazamiento de I2_3
fprintf('\n\nImagen 2_3\nDesplazamiento en filas: %f',dy13)
fprintf('\nDesplazamiento en columnas: %f',dx13)

[dx14 dy14]=LucasKanade(I1,I24,k,4); %Calculamos el desplazamiento de I2_4
fprintf('\n\nImagen 2_4\nDesplazamiento en filas: %f',dy14)
fprintf('\nDesplazamiento en columnas: %f\n',dx14)

close all;

%Muestro Ia=I1(x+da), Ib=I1(x+db), Ra=|I2_1-Ia|, Rb=|I2_1-Ib|
mostrar(I1,I21,dx11,dy11,dx1,dy1,1); 

%Muestro Ia=I1(x+da), Ib=I1(x+db), Ra=|I2_1-Ia|, Rb=|I2_1-Ib|
mostrar(I1,I22,dx12,dy12,dx2,dy2,3);

%Muestro Ia=I1(x+da), Ib=I1(x+db), Ra=|I2_1-Ia|, Rb=|I2_1-Ib|
mostrar(I1,I23,dx13,dy13,dx3,dy3,5);

%Muestro Ia=I1(x+da), Ib=I1(x+db), Ra=|I2_1-Ia|, Rb=|I2_1-Ib|
mostrar(I1,I24,dx14,dy14,dx4,dy4,7);

function mostrar(I,I2,dx11,dy11,dx1,dy1,n)
    a=floor(n/2);
    a=a+1;
    Ia=Interpolar(I,dx1,dy1);
    Ib=Interpolar(I,dx11,dy11);
    Ra=abs(I2-Ia);   Rb=abs(I2-Ib);
    figure(n)
    subplot(1,2,1)
    imagesc(Ib);
    title(['I_b = I_{2\_',num2str(a),'}(x+d_b)']);
    colormap gray;
    subplot(1,2,2)
    imagesc(Ia);
    title(['I_a = I_{2\_',num2str(a),'}(x+d_a)']);
    colormap gray;
    figure(n+1)
    subplot(1,2,2)
    imagesc(Ra);
    title(['R_a = |I_{2\_',num2str(a),'} - I_a|']);
    colormap gray;
    subplot(1,2,1)
    imagesc(Rb);
    title(['R_b = |I_{2\_',num2str(a),'} - I_b|']);
    colormap gray;
end

function [dx dy]=LucasKanade(I1,I2,k,n)
%Programa:      CalcularDespGlobal_Taylor.m
%Descripcion:   Calcular el desplazamiento global entre dos imagenes usando
%               la aproximacion de Taylor de primer orden
%Autor:         FJHL
%Act:           06/Sep/2018
    I1_new=zeros(50,50);
    delta=20; %Para entrar en el while
    acu=0; %Para contar las iteraciones
    dx=0; dy=0; %Distancias a calcular 
    while(delta>0.0001) %Condicion de paro
        acu=acu+1;
        
        %Derivada con respecto del tiempo
        dI2_dt=I2-I1;

        %Derivadas espaciales (Dx, Dy)
        dI1_dx=zeros(50,50);
        dI1_dy=zeros(50,50);
        for x=2:50
            for y=2:50
                dI1_dx(y,x)=I1(y,x)-I1(y,x-1);
                dI1_dy(y,x)=I1(y,x)-I1(y-1,x);
            end
        end

        %Calcular matriz A:
        a11= sum(sum(dI1_dx.^2));
        a12= sum(sum(dI1_dx.*dI1_dy));
        a21= a12;%Ya que la matriz es simetrica
        a22= sum(sum(dI1_dy.^2));
        A=[a11 a12;a21 a22];

        %Calcular el vector de terminos independientes
        b1= -sum(sum(dI1_dx.*dI2_dt));
        b2= -sum(sum(dI1_dy.*dI2_dt));
        b=[b1 b2]';

        %Calcular d=(dx,dy)
        DeltaD=inv(A)*b;

        %Aplicar el desplazamiento encontrado a la imagen I1.
        dx=dx+DeltaD(1); %Calculo las distancias en x
        dy=dy+DeltaD(2); %Calculo las distancias en y
        I1_new=Interpolar(I1,DeltaD(1),DeltaD(2)); %Interpolo la imagen
        delta=sqrt(DeltaD(1)^2+DeltaD(2)^2); %Calulo la magnitud de DeltaD
        I1=I1_new;
        
        if k==1 %Muestro la animacion 
            figure(20)
            subplot(1,2,2)
            imagesc(I1)
            colormap gray;
            title(['I_1(x+d)'])
            subplot(1,2,1)
            imagesc(I2)
            colormap gray;
            title(['I_{2\_',num2str(n)','}(x)'])
            pause(0.1);
        end
    end
end

function I1_new=Interpolar(I1,d1,d2)
    I1_new=zeros(50,50);
    for x=1:50
        for y=1:50
            x_new=x-d1; %La nueva x
            y_new=y-d2; %La nueva y
            if(x_new>=1 && x_new<=50 && y_new>=1 && y_new<=50) %Si se puede mover dentro de la imagen
                I1_new(y,x)=interp2(I1,x_new,y_new); %Interpolamos
            else
                I1_new(y,x)=0;
            end
        end
    end
end

function [dx dy]=FuerzaBruta(I1,I2,n)
    aux=(n-1)/2; 
    desp=-aux; %Desde donde va a iniciar
    for i=1:n %Formo los vectores de desplazamiento
        xd(i)=desp;
        yd(i)=desp;
        desp=desp+1;
    end

    for q=1:n 
        for w=1:n
            acu=0;
            for i=1:50
                for j=1:50
                    if(i-yd(q)>=1&&i-yd(q)<=50&&j-xd(w)>=1&&j-xd(w)<=50) %Si se puede mover dentro de la imagen
                        %Calculo el error para la imagen con un determinado desplazamiento
                        acu=acu+(I2(i,j)-I1(i-yd(q),j-xd(w)))^2; 
                    end
                end
            end
            er(q,w,1)=acu; %Guardo el error de la imagen 
            er(q,w,2)=yd(q); %Guardo su coordenada y
            er(q,w,3)=xd(w); %Guardo su coordenada x
        end
    end
    menor=er(1,1,1); 
    for i=1:n
        for j=1:n
           if(er(i,j,1)<menor) %Si hay otro menor
               menor=er(i,j,1); %Guardo el menor
               dy=er(i,j,2); %Guardo su coordenada y
               dx=er(i,j,3); %Guardo su coordenada x
           end
        end
    end
end