function [win] = slotModel_kCu(win,slotProf,Aslot)

m     = slotProf.m;
q     = slotProf.q;
xMax  = slotProf.xMax;
xMin  = slotProf.xMin;
hSlot = slotProf.hSlot;

kCu      = win.kcu;
ins      = win.condIns;
condType = win.condType;
condHB   = win.condShape;
condHB   = 1;
nCond    = win.nCond;

if strcmp(condType,'Round')
    rMax   = ((Aslot*kCu/nCond)/pi)^0.5;
    wMax   = 2*rMax;
    hMax   = 2*rMax;
    condHB = 1;
else
    wMax = ((Aslot*kCu/nCond)/condHB)^0.5;
    hMax = condHB*wMax;
end


wVector = linspace(wMax*0.1,wMax*1.1,201);
hVector = wVector*condHB;

nVector = zeros(size(wVector));
kVector = zeros(size(wVector));

for ii=1:length(wVector)
    
    indexCond = 0;
    
    hCond = hVector(ii);
    bCond = bVector(ii);
    
    if strcmp(condType,'Round')
        nh = floor(hSlot/(hCond)); % max number of height divisions
    else
        nh = floor(hSlot/(hCond+2*tol)); % max number of height divisions
    end
    
    for xx=1:1:nh
        if xx==1
            x0 = x2-tol-hCond/2;
            xMin = x0-hCond;
            % xMax = x0+hCond;
            yLim = m*xMin+q;
            nw = floor(2*yLim/(wCond+2*tol));
            if rem(nw,2)==0
                flagEven = 1;
            else
                flagEven = 0;
            end
        else
            if strcmp(condType,'Round')
                x0 = x0-2*(rCond+tol)*sin(pi/3);
                xMin = x0-hCond;
                % xMax = x0+hCond;
                yLim = m*xMin+q;
                nw = floor(2*yLim/(wCond+2*tol));
                if flagEven
                    if rem(nw,2)==0
                        nw = nw-1;
                    end
                    flagEven = 0;
                else
                    if rem(nw,2)>0
                        nw = nw-1;
                    end
                    flagEven = 1;
                end
            else
                x0 = x0-hCond-2*tol;
                xMin = x0-hCond;
                % xMax = x0+hCond;
                yLim = m*xMin+q;
                nw = floor(2*yLim/(wCond+2*tol));
                if nw/2==0
                    flagEven = 1;
                else
                    flagEven = 0;
                end
            end
        end
        indexCond=indexCond+nw;
    end
    nVector(ii) = indexCond;
end

nVector(nVector>nCond) = nCond;
nVector(nVector<nCond) = 0;

if strcmp(condType,'Round')
    kVector = nVector.*(wVector/2).^2*pi/Aslot;
else
    kVector = nVector.*wVector.*hVector;
end

[nVector,index] = sort(nVector);
kVector = kVector(index);
wVector = wVector(index);
hVector = hVector(index);

tmp = abs(kVector-kCu);
ii = 1:1:numel(tmp);
index = ii(tmp==min(tmp));


wCond = wVector(index);
hCond = hVector(index);
kCond = kVector(index);
nCond = nVector(index);











