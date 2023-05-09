function [Iw, nBit, MSE] = singleLayerEmbeddingA(I,IX,IY,iPosX,jPosX,iPosY,jPosY,noiseLevel,pFor,Thresh,Payload,Mode)
Iw = I;

[d1,d2] = size(I);

X = unique(noiseLevel);

data = randperm(512^2);
dis = 0;
nBit = 0;
nData = 1;
nData3 = 1;
for iP=67:pFor
    
    
    %Is payload satisfied?
    if(nBit < Payload && noiseLevel(iP)<=Thresh)
        
        index = find(Mode(1,:)==noiseLevel(iP));
        
        switch Mode(2,index)
            case 1
                
                if(abs(IY(iP))==0.5 && abs(IX(iP))==0.5)
                    bit = mod(data(fix(iP/255)*500+mod(iP,255)),3);
                    nData3 = nData3 + 1;
                    nBit = nBit+log2(3);
                    dis = dis + 2/3;
                    if(bit == 1)
                        Iw(iPosX(iP),jPosX(iP)) = I(iPosX(iP),jPosX(iP)) + sign(IX(iP))*1;
                    elseif(bit == 2 )
                        Iw(iPosY(iP),jPosY(iP)) = I(iPosY(iP),jPosY(iP)) + sign(IY(iP))*1;
                    end
                    continue
                end
                
                %     if nData <= length(str2)
                %|x| and |y| = 1.5
                if(abs(IY(iP))==1.5 && abs(IX(iP))==1.5)
                    bit = mod(data(fix(iP/255)*500+mod(iP,255)),2);
                    nData = nData + 1;
                    nBit = nBit+1;
                    dis = dis + 1;
                    if(bit ==1)
                        Iw(iPosX(iP),jPosX(iP)) = I(iPosX(iP),jPosX(iP)) + sign(IX(iP))*1;
                        Iw(iPosY(iP),jPosY(iP)) = I(iPosY(iP),jPosY(iP)) + sign(IY(iP))*1;
                    end
                    continue
                end
                
                %|y| = 0.5
                if(abs(IY(iP))==0.5)
                    bit = mod(data(fix(iP/255)*500+mod(iP,255)),2);
                    nData = nData + 1;
                    nBit = nBit+1;
                    dis = dis + 3/2;
                    if(bit == 0)
                        Iw(iPosX(iP),jPosX(iP)) = I(iPosX(iP),jPosX(iP)) + sign(IX(iP))*1;
                    elseif(bit ==1)
                        Iw(iPosX(iP),jPosX(iP)) = I(iPosX(iP),jPosX(iP)) + sign(IX(iP))*1;
                        Iw(iPosY(iP),jPosY(iP)) = I(iPosY(iP),jPosY(iP)) + sign(IY(iP))*1;
                    end
                    continue
                end
                %
                
                %|x| = 0.5
                if(abs(IX(iP))==0.5)
                    bit = mod(data(fix(iP/255)*500+mod(iP,255)),2);
                    nData = nData + 1;
                    nBit = nBit+1;
                    dis = dis + 3/2;
                    if(bit == 0)
                        Iw(iPosY(iP),jPosY(iP)) = I(iPosY(iP),jPosY(iP)) + sign(IY(iP))*1;
                    elseif(bit ==1)
                        Iw(iPosX(iP),jPosX(iP)) = I(iPosX(iP),jPosX(iP)) + sign(IX(iP))*1;
                        Iw(iPosY(iP),jPosY(iP)) = I(iPosY(iP),jPosY(iP)) + sign(IY(iP))*1;
                    end
                    continue
                end
                
                %     end
                
                dis = dis + 2;
                Iw(iPosX(iP),jPosX(iP)) = I(iPosX(iP),jPosX(iP)) + sign(IX(iP))*1;
                Iw(iPosY(iP),jPosY(iP)) = I(iPosY(iP),jPosY(iP)) + sign(IY(iP))*1;
            case 2
                if(abs(IX(iP))==0.5 )
                    bit = mod(data(fix(iP/255)*500+mod(iP,255)),2);
                    nBit = nBit+1;
                    dis = dis + 0.5;
                    if(bit == 0)
                        Iw(iPosX(iP),jPosX(iP)) = I(iPosX(iP),jPosX(iP)) + sign(IX(iP))*1;
                    elseif(bit == 1 )
                        Iw(iPosX(iP),jPosX(iP)) = I(iPosX(iP),jPosX(iP)) + sign(IX(iP))*1;
                    end
                else
                    dis = dis + 1;
                    Iw(iPosX(iP),jPosX(iP)) = I(iPosX(iP),jPosX(iP)) + sign(IX(iP))*1;
                end
                
                if(abs(IY(iP))==0.5 )
                    bit = mod(data(fix(iP/255)*500+mod(iP,255)),2);
                    nBit = nBit+1;
                    dis = dis + 0.5;
                    if(bit == 0)
                        Iw(iPosY(iP),jPosY(iP)) = I(iPosY(iP),jPosY(iP)) + sign(IY(iP))*1;
                    elseif(bit == 1 )
                        Iw(iPosY(iP),jPosY(iP)) = I(iPosY(iP),jPosY(iP)) + sign(IY(iP))*1;
                    end
                else
                    dis = dis + 1;
                    Iw(iPosY(iP),jPosY(iP)) = I(iPosY(iP),jPosY(iP)) + sign(IY(iP))*1;
                end
                
               
        end
        
    end
    
    %
end

MSE = sum(sum(abs(Iw-I)));
% MSE = dis;
% ps = 10*log10(255^2*d1*d2/MSE);


end