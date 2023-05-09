function [Thresh,bits,MSE,pFor, errorX, errorY, iPosX, jPosX,iPosY,jPosY,noiseLevel2,Mode] = Geosedic2_circleA(I,dir,radius,payload)
[d1,d2] = size(I);

Iper = zeros(1,d1*d2/2);
iPos = zeros(1,d1*d2/2);
jPos = zeros(1,d1*d2/2);
noiseLevel = zeros(1,d1*d2/2);

pfor = 1;
flag = zeros(size(I));
for i=2:d1-1
    if dir+mod(i,2)==2
        k=0;
    else
        k=dir+mod(i,2);
    end
    for j=2+k:2:d2-1
        a = I(i,j-1);
        b = I(i-1,j);
        c = I(i,j+1);
        d = I(i+1,j);
        Iper(pfor) = I(i,j);
        noiseLevel(pfor) = abs(a-c)+abs(b-d)+abs(a-b)+abs(c-d)+abs(a-d)+abs(b-c);
        iPos(pfor) = i;
        jPos(pfor) = j;
        pfor = pfor + 1;
        flag(i,j) = 1;
    end
end
pfor = pfor - 1;
Iper = Iper(1:pfor);
iPos = iPos(1:pfor);
jPos = jPos(1:pfor);
noiseLevel = noiseLevel(1:pfor);

errorX = zeros(1,pfor);
errorY = zeros(1,pfor);

iPosX = zeros(1,pfor);
jPosX = zeros(1,pfor);

iPosY = zeros(1,pfor);
jPosY = zeros(1,pfor);
noiseLevel2 = zeros(1,pfor);

%------------permutation---------
chaos = zeros(1,pfor);
chaos(1)=0.5;
for i=1:pfor-1
    chaos(i+1)=3.7*chaos(i)*(1-chaos(i));
end
[~,ind]=sort(chaos);%将产生的混沌序列进行排序

iPos = iPos(ind);
jPos = jPos(ind);


%---------Geosedic path
pFor = 1;
for i=1:pfor
    X = [];
    xpos = [];
    ypos = [];
    node = 0;
    
    if flag(iPos(i),jPos(i)) == 0
        continue
    end
    
    for ii = iPos(i)-radius:iPos(i)+radius
        for jj = jPos(i)-radius:jPos(i)+radius
            if ii >=1 && ii<=d1 && jj >=1 && jj<=d2
            if flag(ii,jj) == 1 && ii~= iPos(i) && jj~= jPos(i)
                X = [X I(ii,jj)];
                xpos = [xpos ii];
                ypos = [ypos jj];
                node = node + 1;
            end
            end
        end
    end
    X = [I(iPos(i),jPos(i)) X];
    xpos = [iPos(i) xpos];
    ypos = [jPos(i) ypos];
    node = node + 1;
    
    if node < 2
        continue
    end

    ohd = zeros(node,node);
    ohd(1:node,1:node) = Inf;

    for ii = 1:node
        ohd(ii,ii) = 0;
    end
    
    for ii = 1:node
        for jj = 1:node
            distance = sqrt((xpos(ii) - xpos(jj))^2+(ypos(ii) - ypos(jj))^2);
            N1 = [I(xpos(ii)-1,ypos(ii)) I(xpos(ii)+1,ypos(ii)) I(xpos(ii),ypos(ii)-1) I(xpos(ii),ypos(ii)+1)];
            N2 = [I(xpos(jj)-1,ypos(jj)) I(xpos(jj)+1,ypos(jj)) I(xpos(jj),ypos(jj)-1) I(xpos(jj),ypos(jj)+1)];
            Gdelta = sum(abs(N1-N2));
            if ii == jj
                ohd(ii,jj) = 0;
                continue
            end
            ohd(ii,jj) = distance + Gdelta;
        end
    end
    
    mind = dijkstra_distance (node, ohd);
    [~,ind] = sort(mind(1:node));
   
    i3 = xpos(ind(2));
    j3 = ypos(ind(2));
    ii = xpos(ind(1));
    jj = ypos(ind(1));
    
    flag(i3,j3) = 0;
    flag(ii,jj) = 0;
    
      
                      ax = I(ii-1,jj); 
       bx = I(ii,jj-1);                dx = I(ii,jj+1);
                      cx = I(ii+1,jj);

       
                      ay = I(i3-1,j3); 
       by = I(i3,j3-1);                dy = I(i3,j3+1);
                      cy = I(i3+1,j3);

       
     noiseLevel2(pFor) = abs(ax-bx)+abs(bx-cx)+abs(cx-dx)+abs(dx-ax)...
                         + abs(ay-by)+abs(by-cy)+abs(cy-dy)+abs(dy-ay);
    
    errorX(pFor) = I(ii,jj) -  ceil( (I(ii-1,jj) + I(ii,jj-1) + I(ii+1,jj) + I(ii,jj+1))/4);
    errorY(pFor) =  I(i3,j3) - ceil( (I(i3-1,j3) + I(i3,j3-1) + I(i3+1,j3) + I(i3,j3+1))/4);
    iPosX(pFor) = ii;
    jPosX(pFor) = jj; 
    iPosY(pFor) = i3;
    jPosY(pFor) = j3;
    pFor = pFor + 1;
    
end
pFor = pFor - 1;

errorX = errorX(1:pFor);
errorY = errorY(1:pFor);
iPosX = iPosX(1:pFor);
jPosX = jPosX(1:pFor);
iPosY = iPosY(1:pFor);
jPosY = jPosY(1:pFor);
noiseLevel2 = noiseLevel2(1:pFor);
errorX = errorX + 0.5;
errorY = errorY + 0.5;


X = unique(noiseLevel2);
HisEC = zeros(max(X)+1,1);
HisEC_new = zeros(max(X)+1,1);
HisED = zeros(max(X)+1,1);
HisED_new = zeros(max(X)+1,1);
Mode = zeros(2,1500);
EC = [];
ED = [];
n = 0;
%-----------------------
for i=67:pFor
        tt = noiseLevel2(i);
        if(abs(errorX(i))==0.5 && abs(errorY(i))==0.5)
            HisEC_new(tt+1) = HisEC_new(tt+1) + log2(3);
            HisED_new(tt+1) = HisED_new(tt+1) + 2/3;
            HisEC(tt+1) = HisEC(tt+1) + 2;
            HisED(tt+1) = HisED(tt+1) + 1;
            continue
        end
        
        if(abs(errorX(i))==1.5 && abs(errorY(i))==1.5)
            HisEC_new(tt+1) = HisEC_new(tt+1) + 1;
            HisED_new(tt+1) = HisED_new(tt+1) + 1;
            HisED(tt+1) = HisED(tt+1) + 2;
            continue
        end
        
        if(abs(errorX(i))==0.5)
            HisEC_new(tt+1) = HisEC_new(tt+1) + 1;
            HisED_new(tt+1) = HisED_new(tt+1) + 3/2;
            HisEC(tt+1) = HisEC(tt+1) + 1;
            HisED(tt+1) = HisED(tt+1) + 3/2;
            continue
        end
        if(abs(errorY(i))==0.5)
            HisEC_new(tt+1) = HisEC_new(tt+1) + 1;
            HisED_new(tt+1) = HisED_new(tt+1) + 3/2;
            HisEC(tt+1) = HisEC(tt+1) + 1;
            HisED(tt+1) = HisED(tt+1) + 3/2;
            continue
        end
        
            HisED_new(tt+1) = HisED_new(tt+1) + 2;
            HisED(tt+1) = HisED(tt+1) + 2;
    
end

eff = HisED./HisEC;
eff_new = HisED_new./HisEC_new;

%----------------------------
EC = 0;
ED = 0;
n = 1;
for T = X
    if eff_new(T+1) <= eff(T+1)
        EC = EC+HisEC_new(T+1);
        ED = EC+HisED_new(T+1);
        Mode(1,n) = T;
        Mode(2,n) = 1;
    else
        EC = EC+HisEC(T+1);
        ED = EC+HisED(T+1);
        Mode(1,n) = T;
        Mode(2,n) = 2;
    end
    
    if EC >= payload
        break;
    end
    n = n + 1;
end
Mode = Mode(1:2,1:n);
Thresh = T;
bits = EC;
MSE = ED;

end