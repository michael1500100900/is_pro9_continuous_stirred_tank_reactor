clc;
clear all;
close all;

load('dane.mat');
t=dane(:,1);
u=dane(:,2);
y1=dane(:,3);
y2=dane(:,4);
na=2;%Rząd mianownika
nb=2;%Rząd licznika
nk=1;%Opóźnienie
n=length(t);

%*************************************************************************
%---LS estymacja---
phi=zeros(n-nk,na+nb+1); %Inicjalizacja phi stezenia
phi2=zeros(n-nk,na+nb+1); %Inicjalizacja phi temperatury
for i=nk+1:n
    %Przygotowanie macierzy regresji dla y
    for j=1:na
        if(i-j)>0
            phi(i-nk,j)=-y1(i-j);
            phi2(i-nk,j)=-y2(i-j);
        end
    end
    %Przygotowanie macierzy regresji dla u
    for j=0:nb
        if(i-j-nk+1)>0
            phi(i-nk,na+j+1)=u(i-j-nk+1);
            phi2(i-nk,na+j+1)=u(i-j-nk+1);
        end
    end
end
yd1=y1(nk+1:n);
yd2=y2(nk+1:n);

%Estymaty parametrów
theta1= phi\yd1;
theta2=phi2\yd2;
a1=theta1(1:na); %wsp mianownika y1
b1=theta1(na+1:end); %wsp licznika y1
a2=theta2(1:na); %wsp mianownika y2
b2=theta2(na+1:end); %wsp licznika y2

%*************************************************************************
%---Weryfikacja modelu z danymi---
y1_hat=zeros(n,1);
y2_hat=zeros(n,1);
for i=nk+1:n
    y1_hat(i)=0;
    y2_hat(i)=0;
    for j=1:na
        if(i-j)>0
            y1_hat(i)=y1_hat(i)-a1(j)*y1(i-j);
            y2_hat(i)=y2_hat(i)-a2(j)*y2(i-j);
        end
    end
    for j=0:nb
        if(i-j-nk+1)>0
            y1_hat(i)=y1_hat(i)+b1(j+1)*u(i-j-nk+1);
            y2_hat(i)=y2_hat(i)+b2(j+1)*u(i-j-nk+1);
        end
    end
end

%*************************************************************************
%---Wyświetlenie wyników---
fprintf('Oszacowane parametrym wyjścia stężenia metodą LS:\n');
fprintf('Współczynniki mianownika (a):\n');
disp(a1);
fprintf('Współczynniki licznika (b):\n');
disp(b1);
disp('Transmitancja dyskretna z estymat parametrów:\n');
sys1=tf(b1',a1',1)
fprintf('Oszacowane parametry wyjścia temperatury metodą LS:\n');
fprintf('Współczynniki mianownika (a):\n');
disp(a2);
fprintf('Współczynniki licznika (b):\n');
disp(b2);
disp('Transmitancja dyskretna z estymat parametrów:\n');
sys2=tf(b2',a2',1)

%*************************************************************************
%---J_FIT---
%J_fit1=(1-(sum((y1(1:end)-y1_hat(1:end)).^2))/(sum((y1(1:end)-(mean(y1(1:end)))).^2)))*100;
J_fit1=(1-norm(y1(2:end)-y1_hat(2:end))/norm(y1(2:end)-mean(y1(2:end))))*100;%wzor M M Michalka prof PP

%J_fit2=(1-(sum((y2(2:end)-y2_hat(2:end)).^2))/(sum((y2(2:end)-(mean(y2(2:end)))).^2)))*100;
J_fit2=(1-norm(y2(2:end)-y2_hat(2:end))/norm(y2(2:end)-mean(y2(2:end))))*100;%wzor M M Michalka prof PP

%Wyświetlenie wskaźnika J_fit
fprintf('\nWskaźnik J_fit dla stężenia: %.2f%%\n',J_fit1);
fprintf('Wskaźnik J_fit dla temperatury: %.2f%%\n',J_fit2);

%*************************************************************************
%---Wykresy---
figure;
plot(t,y1,'b',t,y1_hat,'r--');
legend('Rzeczywiste','Estymowane');
title('Porównanie stężenia rzeczywistego i estymowanego');
xlim([1 750]);
xlabel('Czas [s]');
ylabel('Stężenie molowe [mol/l]');

figure
plot(t,y2,'b',t,y2_hat,'--');
legend('Rzeczywiste','Estymowane');
title('Porównanie temperatury rzeczywistej i estymowanej');
xlim([1 750]);
xlabel('Czas [s]');
ylabel('Temperatura [K]');