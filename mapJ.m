function [Jx, Jy] = mapJ(JxMask,JyMask,J)
NoE = length(J);
Jx = 0;
Jy = 0;
for i = 1:NoE
    JxIdx = 1;
    JyIdx = 1;
    if(JxMask(i) == 1)
        Jx = [Jx,J(i)];
        JxIdx = JxIdx + 1;
    elseif(JxMask(i) == -1)
        Jx = [Jx,-J(i)];
        JxIdx = JxIdx + 1;
    elseif(JyMask(i) == 1)
        Jy = [Jy,J(i)];
        JyIdx = JyIdx + 1;        
    elseif(JyMask(i) == -1)
        Jy = [Jy,-J(i)];
        JyIdx = JyIdx + 1;
    end
end
Jx = Jx(2:end);
Jy = Jy(2:end);
end