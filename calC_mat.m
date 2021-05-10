function C = calC_mat(indMeta)
% constants
eps0 = 8.8541e-12;
NoI = 20; % number of Image Charges
a = 0; %DNA


% Inductor parameters
OD = indMeta(1);
LW = indMeta(2);
LS = indMeta(3);
N = indMeta(4);
SEP = indMeta(5);
dl = indMeta(6);
epsR = indMeta(7);
eps_sur = 80;

%% map coordinates
[rho_t_x,rho_t_y,NoE] = mapRhoCoord(OD,LW,LS,N,dl);

rho_b_x = -rho_t_x;
rho_b_y = rho_t_y;

% calculate z coordinates of image charge
ImQz = zeros(NoI,1);
ImQz(1) = a;
for idx = 2:NoI
    if mod(idx,2) == 0
        ImQz(idx) = -idx*SEP-a;
    else
        ImQz(idx) = (idx-1)*SEP+a;
    end
end

% calculate image charge factor
ImQfactor = zeros(NoI,1); % 
factor1 = (eps_sur-epsR)/(eps_sur+epsR);
factor2 = (epsR-eps_sur)/(eps_sur+epsR);
for idx = 1:NoI
   if mod(idx,2) == 1
       ImQfactor(idx) = ((-1*factor1*factor2)^((idx-1)/2))*(1-factor1);
   else
       ImQfactor(idx) = ((-factor1)^((idx-2)/2))*((factor2)^(idx/2))*(1-factor1);
   end

%% assign rho
inductorLength = NoE*dl; % total length of inductor
rho0 = 1;
q0 = rho0*inductorLength/pi;
s = linspace(-inductorLength/2,inductorLength/2,NoE);
rho_t = rho0*sin(s*pi/inductorLength); % Top inductor 
rho_b = -rho0*sin(s*pi/inductorLength); % Bottom inductor

%% integrate

invCIntegrand = zeros(NoE,NoE);
invCmIntegrand = zeros(NoE,NoE);

for i = 1:NoE
    for j = 1:NoE
        Rij = calRij(rho_t_x(i),rho_t_y(i),rho_t_x(j),rho_t_y(j));
        Rijtb_rho = calRij(rho_t_x(i),rho_t_y(i),rho_b_x(j),rho_b_y(j));
        Rijtb = sqrt(Rijtb_rho^2 + SEP^2);
        
        % distance btw charge and image charges
        Rij_image = sqrt(Rij^2+(ImQz).^2);
        Rij_tb_image = sqrt(Rijtb_rho^2+(ImQz+SEP).^2);
        if (i == j)
            invCIntegrand(i,j) = rho_t(i)*(1+ImQfactor(1))*rho_t(j)*0.8814*2*(dl+dl);
            invCmIntegrand(i,j) = (rho_t(i)*(1+ImQfactor(1))*rho_b(j)/Rijtb)*(dl^2);
        else
            invCIntegrand(i,j) = (rho_t(i)*(1+ImQfactor(1))*rho_t(j)/(Rij))*(dl^2);
            invCmIntegrand(i,j) = (rho_t(i)*(1+ImQfactor(1))*rho_b(j)/Rijtb)*(dl^2);
        end

        for k = 2:NoI
            invCIntegrand(i,j) = invCIntegrand(i,j) + ...
                (rho_t(i)*ImQfactor(k)*rho_t(j)/Rij_image(k))*(dl^2);
            invCmIntegrand(i,j) = invCmIntegrand(i,j) + ...
                (rho_t(i)*ImQfactor(k)*rho_b(j)/Rij_tb_image(k))*(dl^2);
        end
    end
end

invC_Sum = sum(sum(invCIntegrand));
invCm_Sum = sum(sum(invCmIntegrand));

invC = (1/(4*epsR*pi*eps0*q0^2))*invC_Sum;
invCm = (1/(4*epsR*pi*eps0*q0^2))*invCm_Sum;
C = 2/(invC+invCm);

end