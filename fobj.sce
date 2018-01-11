// Universidade Federal Fluminense
// Departamento de Engenharia Química e de Petróleo
// Codigo para solução do modelo de`Produção Gas Lift
// Alunos: Matheus e Tie
// Prof. Lizando Santos

// Fontes: Alessandro Jacoud Peixoto  Diego Pereira-DiasArthur F. S. Xaud Argimiro Resende Secchi
//        Modelling and Extremum Seeking Control of Gas Lifted Oil Wells
//        IFAC-PapersOnLine 48-6 (2015) 021–026


//        Model Predictive Control with quality requirements on petroleum production platforms
//        C.H.P. Ribeiro, S.C. Miyoshi, A.R. Secchi, A. Bhaya

function [fobj] = fobj(par)



function [f]=fun(t,x,par)  // Modelo Cido e Peixoto (Original - sem per)

// Parametros

M = 0.0289; // Massa molar do gas [kg/mol]
Ta = 293; // Temperatura da regiao anular [K]
La	= 230.87; // Comprimento da regiao anular [m]
Va	= 29.012; // Volume da regiao anular [m3]
ro	= 923.9;  // Densidade do óleo no reservatorio [kg/m3]
Tw	= 293;  // Temperatura no tubo [K]
Lw	= 1217;  // Comprimento do tubo [m]
Lr	= 132; // Distancia do reservatorio ate o ponto de injeção [m]
Ar	= 0.203 ; // Area da seção transversal abaixo do ponto injecao [m2]
Aw  = 0.203; // Area da seção transversal acima do ponto injecao [m2]
Cpc	= 1.655e-3; // [m2]
Civ	= 15e-5; // [m2]
Cr = 2.623e-4; // [m2]
Pr = 2.5497295e7// Pressao no reservatorio longe da cabeça do poço [Pa];
Ps = 3.704669e6 // Pressao no manifold [Pa];
vo = 1/ro;  // volume específico do óleo
g = 9.81; // m/s2
R = 8.314; // J/mol.K

// Pressoes

Pai=((R*Ta/(Va*M)+g*La/Va)*x(1)); // Pressao na regiao anular (na injeção)
Pwh=((R*Tw/M)*(x(2)/(Lw*Aw-vo*x(3) ))); // Pressao na cabeça do poço
Pwi=(Pwh + (g/Aw)*(x(2)+x(3) ) - (ro*Lr*Ar)); // Pressao da coluna no ponto de injeção
Pwb=(Pwi + ro*g*Lr); // Pressao no fundo do poço

// Densidades

ro_ai = (M*Pai/(R*Ta) + g*La/Va)*x(1); //Densidade do gás na regiao anular
ro_m = (( (x(2) + x(3) )  - ro*Lr*Ar )/ (Lw*Aw) ); // densidade na cabeça do poço

// Vazoes

y1 = max(0,(Pai-Pwi));
y2 = max(0,(Pr-Pwb));
y3 = max(0,(Pwh-Ps));

wiv=(Civ*sqrt(ro_ai*y1));  //vazão de gás que sai da regiao anular para o tubo
wro = (Cr*sqrt(ro*y2)); // vazão de gás do reservatorio para o tubo
wpc = Cpc*sqrt(ro_m*y3)  // vazao na cabeça do poço
wpg = x(2)*wpc/( x(2) + x(3) ); // vazão de gás na cabeça do poço
wpo = x(3)*wpc/( x(2) + x(3) ); // vazão de óleo na cabeça do poço

u = par; // Vazao de gas de injecao Kg/s

// Equacoes Diferenciais

dx(1)=(u - wiv) ;
dx(2)=(wiv + wro*0.0 - wpg) ;
dx(3)=(wro - wpo) ;

f = [dx];
endfunction



function [f]=modelo(x)  // Modelo (Esta função é para PLOTAR as variaveis)

// Parametros

M = 0.0289; // Massa molar do gas [kg/mol]
Ta = 293; // Temperatura da regiao anular [K]
La	= 230.87; // Comprimento da regiao anular [m]
Va	= 29.012; // Volume da regiao anular [m3]
ro	= 923.9;  // Densidade do óleo no reservatorio [kg/m3]
Tw	= 293;  // Temperatura no tubo [K]
Lw	= 1217;  // Comprimento do tubo [m]
Lr	= 132; // Distancia do reservatorio ate o ponto de injeção [m]
Ar	= 0.203 ; // Area da seção transversal abaixo do ponto injecao [m2]
Aw  = 0.203; // Area da seção transversal acima do ponto injecao [m2]
Cpc	= 1.655e-3; // [m2]
Civ	= 15e-5; // [m2]
Cr = 2.623e-4; // [m2]
Pr = 2.5497295e7// Pressao no reservatorio longe da cabeça do poço [Pa];
Ps = 3.704669e6 // Pressao no manifold [Pa];
vo = 1/ro;  // volume específico do óleo
g = 9.81; // m/s2
R = 8.314; // J/mol.K

// Pressoes

Pai=((R*Ta/(Va*M)+g*La/Va)*x(1)); // Pressao na regiao anular (na injeção)
Pwh=((R*Tw/M)*(x(2)/(Lw*Aw-vo*x(3) ))); // Pressao na cabeça do poço
Pwi=(Pwh + (g/Aw)*(x(2)+x(3) ) - (ro*Lr*Ar)); // Pressao da coluna no ponto de injeção
Pwb=(Pwi + ro*g*Lr); // Pressao no fundo do poço

// Densidades

ro_ai = (M*Pai/(R*Ta) + g*La/Va)*x(1); //Densidade do gás na regiao anular
ro_m = (( (x(2) + x(3) )  - ro*Lr*Ar )/ (Lw*Aw) ); // densidade na cabeça do poço

// Vazoes

y1 = max(0,(Pai-Pwi));
y2 = max(0,(Pr-Pwb));
y3 = max(0,(Pwh-Ps));

wiv=(Civ*sqrt(ro_ai*y1));  //vazão de gás que sai da regiao anular para o tubo
wro = (Cr*sqrt(ro*y2)); // vazão de gás do reservatorio para o tubo
Pwb 

wpc = Cpc*sqrt(ro_m*y3)  // vazao na cabeça do poço
wpg = x(2)*wpc/( x(2) + x(3) ); // vazão de gás na cabeça do poço
wpo = x(3)*wpc/( x(2) + x(3) ); // vazão de óleo na cabeça do poço

u = par; // Vazao de gas de injecao Kg/s

f = [wiv; wro; wpc; wpg; wpo; Pai; Pwh; Pwi; Pwb ];

endfunction

k = 1;
//for j =3:0.5:7
//
//par=j;  // Vazao de gas de injecao Kg/s
//u(k) = par
t0 = 0;  // Instante inicial [s]
//tf = 3600*10;  // Instante final [s]
tf = 5;  // Instante final [s]
t = [t0:1:tf]; // Intervalo de tempo [s]
//y0= [4350.1;10951;86038];
y0 = [4350.1;10951;86038]; // Condicçoes iniciais da massa [kg]

x =ode(y0,t0,t,fun); // Integração dos sistema de EDO

// Funcao para plotar as variaveis (chamo outra funcao com o modelo)
dimension = size(x);
dim = dimension(2);
for i=1:dim
    vazoes(:,i) = modelo(x(:,i));  // calcula a vazao para cada instante de tempo
end

// Funções para plotar


//figure(1)
//subplot(311),plot2d(t,vazoes(1,:));  //vazão de gás que sai da regiao anular para o tubo
//xtitle("Vazão de gás regiao anular - [kg/s] x tempo [s]")
//subplot(312),plot2d(t,vazoes(2,:));  // vazão de gás do reservatorio para o tubo
// xtitle("Vazão de gás do reservatorio para o tubo - [kg/s] x tempo [s]")
//subplot(313),plot2d(t,vazoes(3,:));  // vazao na cabeça do poço
//xtitle("Vazão total na cabeça do poço - [kg/s] x tempo [s]")
 
 
//figure(2)
//subplot(311),plot2d(t,vazoes(4,:));  // vazão de gás na cabeça do poço
//xtitle("Vazão de gás na cabeça do poço - [kg/s] x tempo [s]")
//subplot(312),plot2d(t,vazoes(5,:));  // vazão de óleo na cabeça do poço
//xtitle("Vazão de óleo na cabeça do poço - [kg/s] x tempo [s]")
//subplot(313),plot2d(t,vazoes(6,:));  // Pressao na regiao anular (na injeção)
//xtitle("Pressão na regiao anular (na injeção) - [Pa] x tempo [s]")
//
//figure(3)
//subplot(311),plot2d(t,vazoes(7,:));  // Pressao na cabeça do poço
//xtitle("Pressão na cabeça do poço - [Pa] x tempo [s]")
//subplot(312),plot2d(t,vazoes(8,:));  // Pressao da coluna no ponto de injeção
//xtitle("Pressão da coluna no ponto de injeção - [Pa] x tempo [s]")
//subplot(313),plot2d(t,vazoes(9,:));  // Pressao no fundo do poco
//xtitle("Pressão no fundo do poço - [Pa] x tempo [s]")
kk = length(vazoes(9,:))
//
////dP(k) = vazoes(5,kk) - 3.704669e6;
//romedio(k) = vazoes(5,kk)
////efect(k) = dP(j)*romedio(j)
////ratio(k) = x(3,kk)/(x(2,kk)+x(3,kk))
//
//k = k+1
//
//end
//figure(1)
//plot(u,romedio)
//figure(2)
//plot(u,ratio)

fobj = -vazoes(5,kk)

endfunction
//
