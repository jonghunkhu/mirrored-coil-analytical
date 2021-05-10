% this function maps the coordinates for a rectangular spiral inductor
function [x,y,NoE] = mapRhoCoord(OD,LW,LS,N,dl)
x = 0;
y = 0;
for nIdx = 1:N
    % Top left corner
    TLcornerXcoord = -OD/2+(nIdx-1)*(LW+LS); 
    TLcornerYcoord = OD/2-(nIdx-1)*(LW+LS); 
    
    X1 = 0;
    Y1 = TLcornerYcoord-LW/2;
    
    X2 = -TLcornerXcoord-LW/2;
    Y2 = TLcornerYcoord-LW/2;
    
    X3 = -TLcornerXcoord-LW/2;
    Y3 = -TLcornerYcoord+LW/2;
    
    X4 = TLcornerXcoord+LW/2;
    Y4 = -TLcornerYcoord+LW/2;
    
    X5 = TLcornerXcoord+LW/2; 
    Y5 = TLcornerYcoord-3*LW/2-LS;
    
    X6 = 0-dl;
    Y6 = TLcornerYcoord-3*LW/2-LS;
    
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
    
    x = [x, n12coordX,n23coordX,n34coordX,n45coordX,n5endcoordX];
    y = [y, n12coordY,n23coordY,n34coordY,n45coordY,n5endcoordY];
end

%{
    if (nIdx == N)
        x = [x, n12coordX,n23coordX,n34coordX,n45coordX,X5];
        y = [y, n12coordY,n23coordY,n34coordY,n45coordY,Y5];
    else
        x = [x, n12coordX,n23coordX,n34coordX,n45coordX,n5endcoordX];
        y = [y, n12coordY,n23coordY,n34coordY,n45coordY,n5endcoordY];
    end
end
%}
x = x(2:end);
y = y(2:end);
NoE = numel(x);
end