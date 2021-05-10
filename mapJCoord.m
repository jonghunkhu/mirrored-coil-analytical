% this function maps the coordinates for a rectangular spiral inductor
function [Jx_xCoord,Jx_yCoord,Jy_xCoord,Jy_yCoord,Jxmask,Jymask,NoE] = mapJCoord(OD,LW,LS,N,dl)

Jx_xCoord = 0;
Jx_yCoord = 0;
Jy_xCoord = 0;
Jy_yCoord = 0;
Jxmask = 0;
Jymask = 0;

for nIdx = 1:N
    startXcoord = -OD/2+(nIdx-1)*(LW+LS);
    startYcoord = OD/2-(nIdx-1)*(LW+LS);
    
    % X1 = startXcoord+LW/2;
    X1 = 0;
    Y1 = startYcoord-LW/2;
    
    X2 = -startXcoord-LW/2;
    Y2 = startYcoord-LW/2;
    
    X3 = -startXcoord-LW/2;
    Y3 = -startYcoord+LW/2;
    
    X4 = startXcoord+LW/2;
    Y4 = -startYcoord+LW/2;
    
    X5 = startXcoord+LW/2;
    Y5 = startYcoord-3*LW/2-LS;
    
    %X6 = -OD/2+(nIdx)*(LW+LS)+LW/2-dl;
    X6 = 0 - dl;
    Y6 = startYcoord-3*LW/2-LS;
    
    n12 = round(abs(X2-X1)/dl);
    n23 = round(abs(Y3-Y2)/dl);
    n34 = round(abs(X4-X3)/dl);
    n45 = round(abs(Y5-Y4)/dl);
    n5end = round(abs(X6-X5)/dl);
    
    n12coordX = linspace(X1,X2-dl,n12);
    n23coordX = linspace(X2,X3,n23);
    n34coordX = linspace(X3,X4+dl,n34);
    n45coordX = linspace(X4,X5,n45);
    n5endcoordX = linspace(X5,X6,n5end);
    
    n12coordY = linspace(Y1,Y2,n12);
    n23coordY = linspace(Y2,Y3+dl,n23);
    n34coordY = linspace(Y3,Y4,n34);
    n45coordY = linspace(Y4,Y5-dl,n45);
    n5endcoordY = linspace(Y5,Y6,n5end);
    
    n12ones = ones(1,n12);
    n23ones = ones(1,n23);
    n34ones = ones(1,n34);
    n45ones = ones(1,n45);
    n5ones = ones(1,n5end);
    
    n12zeros = zeros(1,n12);
    n23zeros = zeros(1,n23);
    n34zeros = zeros(1,n34);
    n45zeros = zeros(1,n45);
    n5zeros = zeros(1,n5end);
    
    % Assign coordinates
%     if (nIdx == N)
%         Jx_xCoord = [Jx_xCoord,n12coordX,n34coordX,X5];
%         Jx_yCoord = [Jx_yCoord,n12coordY,n34coordY,Y5];
%         Jxmask = [Jxmask,n12ones,n23zeros,-1*n34ones,n45zeros,1];
%         
%         Jy_xCoord = [Jy_xCoord,n23coordX,n45coordX];
%         Jy_yCoord = [Jy_yCoord,n23coordY,n45coordY];
%         Jymask = [Jymask,n12zeros,-1*n23ones,n34zeros,n45ones,0];
    
        Jx_xCoord = [Jx_xCoord,n12coordX,n34coordX,n5endcoordX];
        Jx_yCoord = [Jx_yCoord,n12coordY,n34coordY,n5endcoordY];
        Jxmask = [Jxmask,n12ones,n23zeros,-1*n34ones,n45zeros,n5ones];
        
        Jy_xCoord = [Jy_xCoord,n23coordX,n45coordX];
        Jy_yCoord = [Jy_yCoord,n23coordY,n45coordY];
        Jymask = [Jymask,n12zeros,-1*n23ones,n34zeros,n45ones,n5zeros];
    
end
Jx_xCoord = Jx_xCoord(2:end);
Jx_yCoord = Jx_yCoord(2:end);
Jy_xCoord = Jy_xCoord(2:end);
Jy_yCoord = Jy_yCoord(2:end);
Jxmask = Jxmask(2:end);
Jymask = Jymask(2:end);
NoE = numel(Jx_xCoord)+numel(Jy_xCoord);
end