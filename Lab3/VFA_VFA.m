%% Este scrip sirve para graficar las respuestas en frecuencia y diagramas
%  de Bode del circuito VFA-VFA a modo comparativo con lo obtenido en
%  LTspice

clear all; close all; clc;

%Amplificador
%%
Ad0 = 100000;                               %Ganancia
f1=10;                                      %Polos
f2=5.06e6; 
w1=2*pi*f1;
w2=2*pi*f2; 
z=[];                                       %Ceros
pol=[-w1 -w2];
gan=Ad0*w1*w2;
[numopamp , denopamp] = zp2tf(z,pol,gan);
tfopamp=tf(numopamp,denopamp);              %Funcion de transferencia

%Grafico amplificador
%%
figure (1);
hold on;
h = bodeplot(tfopamp);
p = getoptions(h);
p.PhaseVisible='on';
p.FreqUnits='Hz';
p.Grid='on';
p.Xlim=[0.1 100000000];
setoptions(h,p);


%%
%{
Ampificador compuesto con segundo AO considerado ideal:
Caclulo de ganancia ideal del segundo opamp considerado como ideal
para lograr un margen de fase de 65 grados.
%}
syms wg;
M_phi_ideal = ((65.0*pi)/180) == pi - atan(wg/w1) - atan(wg/w2);
wg_ideal = solve(M_phi_ideal);                              %Resuelve
wg_ideal = double(wg_ideal);                                %Lo coloca con decimales
fg_ideal = wg_ideal/(2*pi);                                 %Encuentra la frecuencia fg ideal

Avfglobalideal=10;                                          %Avfglobalideal=1/kglobal
Avf2sinpol=((Avfglobalideal)*wg_ideal)/(Ad0*w1)             %Avf2 resultante de opamp2 no tiene polos

%%
%{
Ampificador compuesto con segundo opamp considerado ideal:
Graficos de:
   Ganancia de lazo abierto del compuesto, 
   Ganacia de lazo cerrado global ideal del compuesto
   Gancia Avf2 considerando al segundo amplificador operacional ideal
%}

Avcomp0=Ad0*Avf2sinpol;                     %Ganacia a lazo abierto en DC del Compuesto (veces) -> (127dB)
z=[];                                       %zeros a lazo abierto del compuesto
pol=[-w1 -w2];                              %polos a lazo abierto del compuesto
gan=Avcomp0*w1*w2;                          %expresion de ganacia compatible con matlab
[numavcomp , denavcomp] = zp2tf(z,pol,gan);
tfavcomp=tf(numavcomp,denavcomp);
h = bodeplot(tfavcomp);                     
p = getoptions(h);
p.PhaseVisible='on';
p.FreqUnits='Hz';
p.Grid='on';
p.Xlim=[0.1 100000000];
setoptions(h,p);


tfAvf2sinpol=tf(Avf2sinpol);
bode(tfAvf2sinpol);                        
tfAvfglobalideal=tf(Avfglobalideal);
bode(tfAvfglobalideal);                    

%%
%Funcion de Transferencia de lazo cerrado del compuesto
kglobal=1/Avfglobalideal;
numAvf=[Ad0*Avf2sinpol*w1*w2];
denAvf=[1 (w1+w2)  (w1*w2+kglobal*Ad0*Avf2sinpol*w1*w2)];
tfAvfcomp=tf(numAvf,denAvf);
bode(tfAvfcomp);
[sos,g] = tf2sos(numAvf,denAvf); %Descomponemos en bicuadrativas
wp=sqrt(sos(1,6));
fp=wp/2*pi;
%qp debe ser aproximadamente 0.707 (1/sqrt(2)) para margen de fase de 65
%grados
qp=1/(sos(1,5)/wp);
%relacion entre wg/wp debe ser arpoximadamente 0.65 para margen de fase de
%65 grdos
r1=wg_ideal/wp;

%Calculo del wh suponiendo el segundo amplificador sin polos
syms wh
BW_ideal = 2 == (sos(1,6)-wh^2)^2+sos(1,5)^2;
wh_ideal = solve(BW_ideal);
wh_ideal = double(wh_ideal);
wh_id= abs(wh_ideal(1));
%Relacion entre wp y ancho de banda wh debe ser aproximadamente 1 para
%margen de fase de 65 gradis
r2=wh_id/wp;
%%

