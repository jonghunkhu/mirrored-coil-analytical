function L = calL(indMeta)
% constants
mu0 = 4*pi*1e-7;
% Inductor parameters
OD = indMeta(1);
LW = indMeta(2);
LS = indMeta(3);
N = indMeta(4);
SEP = indMeta(5);
dl = indMeta(6);
%% map coordinates
[Jx_t_x,Jx_t_y,Jy_t_x,Jy_t_y,Jxmask,Jymask,NoE] = mapJCoord(OD,LW,LS,N,dl);
Jx_b_x = -Jx_t_x;
Jx_b_y = Jx_t_y;
Jy_b_x = -Jy_t_x;
Jy_b_y = Jy_t_y;
xNoE = length(Jx_t_x);
yNoE = length(Jy_t_x);

% figure;
% hold on;
% scatter(Jx_t_x,Jx_t_y);
% scatter(Jy_t_x,Jy_t_y);
% 
% figure;
% hold on;
% scatter(Jx_b_x,Jx_b_y);
% scatter(Jy_b_x,Jy_b_y);

%% assign j
inductorLength = NoE*dl; % total length of inductor
s = linspace(-inductorLength/2,inductorLength/2,NoE);
I0 = 1;
I = I0*cos(s*pi/inductorLength);
J = I; % current to current density
[Jx_t,Jy_t] = mapJ(Jxmask,Jymask,J);
Jx_b = Jx_t;
Jy_b = -Jy_t;
%%
JxIntegrand = zeros(xNoE,xNoE);
JxmIntegrand = zeros(xNoE,xNoE);
for i = 1:xNoE
    for j = 1:xNoE
        Rij = calRij(Jx_t_x(i),Jx_t_y(i),Jx_t_x(j),Jx_t_y(j));
        Rijtb_rho = calRij(Jx_t_x(i),Jx_t_y(i),Jx_b_x(j),Jx_b_y(j));
        Rijtb = sqrt(Rijtb_rho^2 + SEP^2);
        if (i == j)
            JxIntegrand(i,j) = Jx_t(i)*Jx_t(j)*0.8814*2*(dl+dl);
        else
            JxIntegrand(i,j) = Jx_t(i)*Jx_t(j)/(Rij)*(dl^2); 
        end
        JxmIntegrand(i,j) = Jx_t(i)*Jx_b(j)/Rijtb*(dl^2);
    end
end

JyIntegrand = zeros(yNoE,yNoE);
JymIntegrand = zeros(yNoE,yNoE);
for i = 1:yNoE
    for j = 1:yNoE
        Rij = calRij(Jy_t_x(i),Jy_t_y(i),Jy_t_x(j),Jy_t_y(j));
        Rijtb_rho = calRij(Jy_t_x(i),Jy_t_y(i),Jy_b_x(j),Jy_b_y(j));
        Rijtb = sqrt(Rijtb_rho^2 + SEP^2);
        if (i == j)
            JyIntegrand(i,j) = Jy_t(i)*Jy_t(j)*0.8814*2*(dl+dl);
        else
            JyIntegrand(i,j) = (Jy_t(i)*Jy_t(j)/(Rij))*(dl^2);
        end
        JymIntegrand(i,j) = (Jy_t(i)*Jy_b(j)/Rijtb)*(dl^2);
    end
end

JxSum = sum(sum(JxIntegrand));
JxmSum = sum(sum(JxmIntegrand));
JySum = sum(sum(JyIntegrand));
JymSum = sum(sum(JymIntegrand));

L = 2*mu0/(4*pi*I0^2)*(JxSum+JySum+JxmSum+JymSum);
end