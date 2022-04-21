%GABRIEL FERNANDES ROCATELLI 
%IMPLEMENTAÇÃO MÉTODO ELEMENTOS FINITOS

%% DADOS
clc;clear all;
elementnumber = 4;
nos = elementnumber + 1;
k1 = 2; %N/mm
k2 = 4; %N/mm
L = 1;
Le = 1/4;
E = [2 2 2 2 2]; %MPa
A = [5 5 5 5]; %mm2
mrigidez = [1 -1;-1 1];
%inc = [0 1/4; 1/4 1/2; 1/2 3/4;3/4 1];

%% INCIDENCIA NODAL 
for i=1:elementnumber
inci(i,:)=[ i  i+1];    %identificação da localização de cada elemento - elemento 1 nós 1 e 2; elemento 2 nós 2 e 3; elemento 3 nós 3 e 4; elemento 4 nós 4 e 5.
end

%% METODO DIFERENCIAL
syms csi 
Nii = [(1-csi)/2];
Njj = [(1+csi)/2]; 
Ni = [diff(Nii)];%ao inves de usar o diferencial, é mais interessante usar a evidencia de valores dentro da matriz calculada e implementar os termos simplificados
Nj = [diff(Njj)];



%% MATRIZ LOCAL
for i=1:elementnumber
    %Ke(:,:,i) = int(2*A(i)*E(i)*(2/Le)*Ni(i)*Nj(i),-1,1);
    K(:,:,i) = (2*A(i)*E(i)*(2/L))*mrigidez;
    if i==1 ;
        K(1,1,1)=K(1,1)+k1; %condição de contorno na extremidade esquerda da barra
    end
    if i==4;
        K(2,2,4)=K(2,2)+k2; %condição de contorno na extremidade direita da barra
    end
end
%% MATRIZ GLOBAL
KG = zeros(nos,nos); %matriz zero que representa a matriz global em dimensões
 for e=1:elementnumber
     loc=[inci(e,1) inci(e,2)]; %posição na qual os elementos entrarão na matriz de rigidez global
     KG(loc,loc) = KG(loc,loc) + K(:,:,e); %matriz global
 end
%% CARREGAMENTO
k=1/4; %coordenada x2 do nó
j=0; %coordenada x1 do nó
for i=1:elementnumber
    T(i)=(10*((j/2) - ((j*csi)/2) + (k/2) + ((k*csi)/2))) - 2; 
    k=k+1/4;
    j=j+1/4;
end
%uma vez com o T avaliado nos índices dos elementos, pode-se implementar a
%força aplicada nos nós a partir da aproximação de Galerkin

%%
%APLICAÇÃO DAS FORÇAS
for i=1:elementnumber
    Fi(i)=(int(T(i)*Nii*(1/8),-1,1)); %laço que calcula a força aplicada na posição Ni.
    Fj(i)=(int(T(i)*Njj*(1/8),-1,1)); %laço que calcula a força aplicada na posição Nj.
    if i==4
        Fj(i)=Fj(i) - k2 * 0.1 * 1; %condição de contorno na extremidade direita
    end
end

%%
%EQUILÍBRIO DE FORÇAS APLICADAS NOS NÓS
F=zeros(5,1); %matriz que representa as dimensões da força em cada nó. Implementada como matriz ao invés de vetor
F(1,1)=Fi(1);
F(2,1)=Fi(2) + Fj(1);
F(3,1)=Fi(3) + Fj(2);
F(4,1)=Fi(4) + Fj(3);
F(5,1)= Fj(4);
%% 
%POR FIM, OS DESLOCAMENTOS PODEM SER LOCALIZADOS ATRAVÉS DO SISTEMA LINEAR:

format long
U = KG\F;
U;

% %FUNÇÕES NODAIS
 syms x
 
 for i=1:elementnumber
     N(i,:) = [i-(x/Le),4*x+1-i]; %declaração da função N(x)
 end
  syms F



% Deformação no elemento para um dado valor de x
for e=1:elementnumber
loc=[inci(e,1) inci(e,2)];
U = U(inci);
U_element = U(loc);
U_element = U_element';
N_element = [N(e,1) N(e,2)];
u(e,:) = N_element*U_element;
coord = [0 0.25 0.5 0.75 1];
end

u_plot_1 = [269/2720 43/102];
x = coord;
y_plot_1 = polyval(u_plot_1,x);
y_plot_1 = y_plot_1';

u_plot_2 = [47/544 6931/16320];
y_plot_2 = polyval(u_plot_2,x);
y_plot_2 = y_plot_2';

u_plot_3 = [31/2720 7543/16320];
y_plot_3 = polyval(u_plot_3,x);
y_plot_3 = y_plot_3';

u_plot_4 = [17759/32640 133/1360];
y_plot_4 = polyval(u_plot_4,x);
y_plot_4 = y_plot_4';


U_plot = zeros(20,1)

for i=1:5
    U_plot(i) = U_plot(i) + y_plot_1(i)
end

for i=6:10
    U_plot(i) = U_plot(i) + y_plot_2(i-5)
end
for i=11:15
    U_plot(i) = U_plot(i) + y_plot_3(i-10)
end
for i=16:20
    U_plot(i) = U_plot(i) + y_plot_4(i-15)
end

figure
plot(U_plot)

% [F]= gauss_quadra(@(x) u,-1,1,6)
% %%Tensão do elemento para um dado valor de x
% for i=1:elementnumber
% 
%     Tensao(i) = diff(u(i)); 
%     Tensao(i) = Tensao(i)*E(i);
% 
% end
% coord=[0 0.25 0.50 0.75 1.0;]
% U = KG\F;
% figure
% plot(coord, U)
% plot(coord, Tensao)
% %professor roberto  daledone machado marcos ardnt 

% x1=0:1;        
% for i=1:length(x1)
%     F=subs(F,x1(i))
% end
