clear all
clc

imgfile = ['J:\Matlab_Program\MATLAB\My_Histogram\PVO_based\data\'];
imgdir = dir([imgfile,'\*.bmp']);
fid=fopen('fileName.txt','wt');
performance = zeros(length(imgdir)*2,100);
location_map = zeros(length(imgdir),100);
tic
for i_img = 5:5

    i_img

%     if i_img == 5 || i_img == 8
%         continue
%     end
    
    if i_img == 2
        stepSize = 1000;
    else
        stepSize = 2000;
    end
    
img = 2*(i_img-1)+1;
I = double(imread([imgfile,'\',imgdir(i_img).name]));

[d1,d2] = size(I);
nIndex = 1;

for payload = [10000]
    bestMode = [];

    Iw = I;  
    dir = 0;
    embedBits = 0;
    dis = 0;
    Re = zeros(50,4);
    tmp = 0;
    maxPSNR = 0;
    for radius = 2:2
        [bin_LM bin_LM_len I] = LocationMap(I);
        [Thresh,bits,MSE,pFor,errorX, errorY, iPosX, jPosX,iPosY,jPosY,noiseLevel,Mode] = Geosedic2_crossA(I,dir,radius,payload*0.5+bin_LM_len);
        if bits >= payload*0.5
            tmp = tmp + 1;
            Re(tmp,1) = bits;
            Re(tmp,2) = 10*log10(255^2*512^2/MSE);
            Re(tmp,3) = Thresh;
            Re(tmp,4) = radius;
            if Re(tmp,2) > maxPSNR
                maxPSNR = Re(tmp,2);
                iPosX2 = iPosX;
                jPosX2 = jPosX;
                iPosY2 = iPosY;
                jPosY2 = jPosY;
                pFor2 = pFor;
                IX = errorX;
                IY = errorY;
                noiseLevel2 = noiseLevel;
                bestMode = Mode;
            end
        end
        
    end
    [val,ind] = max(Re(:,2));
    Thresh = Re(ind,3);
    radius = Re(ind,4);
    if Re(ind,1) < payload*0.5
        break
    end
    [Iw, nBit, MSE] = singleLayerEmbeddingA(I,IX,IY,iPosX2,jPosX2,iPosY2,jPosY2,noiseLevel2,pFor2,Thresh,payload*0.5,bestMode);
    dis = dis + MSE;
    embedBits = embedBits+ nBit;
    
    %---------2nd layler
    
    dir = 1;

    Re2 = zeros(50,4);
    tmp = 0;
    maxPSNR = 0;
    for radius = 2:2
        [bin_LM bin_LM_len Iw] = LocationMap_circle(Iw);
        [Thresh,bits,MSE,pFor, errorX, errorY, iPosX, jPosX,iPosY, jPosY,noiseLevel,Mode] = Geosedic2_circleA(Iw,dir,radius,payload*0.5+bin_LM_len);
        if bits >= payload*0.5
            tmp = tmp + 1;
            Re2(tmp,1) = bits;
            Re2(tmp,2) = 10*log10(255^2*512^2/MSE);
            Re2(tmp,3) = Thresh;
            Re2(tmp,4) = radius;
            if Re2(tmp,2) > maxPSNR
                maxPSNR = Re2(tmp,2);
                iPosX2 = iPosX;
                jPosX2 = jPosX;
                iPosY2 = iPosY;
                jPosY2 = jPosY;
                pFor2 = pFor;
                IX = errorX;
                IY = errorY;
                noiseLevel2 = noiseLevel;
                bestMode = Mode;
            end
        end
    end
        
    [val,ind] = max(Re2(:,2));
    Thresh = Re2(ind,3);
    radius = Re2(ind,4);
    if Re2(ind,1) < payload*0.5
        break
    end
    [Iw, nBit, MSE] = singleLayerEmbedding2A(Iw,IX,IY,iPosX2,jPosX2,iPosY2,jPosY2,noiseLevel2,pFor2,Thresh,payload*0.5,bestMode);
    dis = dis + MSE;
    embedBits = embedBits+ nBit;
    dis2 = sum(sum((I-Iw).^2));
    if embedBits < payload
        break
    end
    performance(img,nIndex) = embedBits;
    performance(img+1,nIndex) = 10*log10(255^2*512^2/dis);
%     if performance(img,nIndex) < payload
%     performance(img,nIndex) = 0;
%     performance(img+1,nIndex) = 0;
%     break
%     end
    nIndex = nIndex + 1;
    
end
end
toc
% save M2.mat performance
% 

