% Pengqunnie
clear; 
close all;
clc;
%% %%%%%%%%%% Airlight Estimation %%%%%%%%%%% %%
I0 = imread('Ir02_0.jpg');
I90 = imread('Ir02_90.jpg');

showImg = figure, imshow(I0,[]);
title('left click to pick up the most polarized pixel in infinite region ')
[x, y] = ginput(1);
inf0 = I0(round(y), round(x), :);
inf0 = double(inf0) -1;
inf0 = min(inf0, 255);

inf90 = I90(round(y), round(x), :);
inf90 = double(inf90) -1;
inf90 = min(inf90, 255);
close(showImg);

A_infi = inf0 + inf90;
DoP_infi = (inf90 - inf0)./A_infi;

%% %%%%%%%%%% Transimission & Dehazing%%%%%%%%%% %%
I1 = double(I0);
I2 = double(I90);
I_tot = I1 + I2;
DoP = (I2 - I1)./I_tot;
DoP = DoP(1,1);
DoP_infi = DoP_infi(1,1);
A_infi = A_infi(1,1);

[m,n,t] = size(I_tot);
Trans(:,:,1)=ones(m,n)-I_tot(:,:,1)*DoP./(A_infi*DoP_infi);
Trans(:,:,2)=ones(m,n)-I_tot(:,:,2)*DoP./(A_infi*DoP_infi);
Trans(:,:,3)=ones(m,n)-I_tot(:,:,3)*DoP./(A_infi*DoP_infi);

L_obj(:,:,1)=I_tot(:,:,1)*(1-DoP/DoP_infi)./Trans(:,:,1);
L_obj(:,:,2)=I_tot(:,:,2)*(1-DoP/DoP_infi)./Trans(:,:,2);
L_obj(:,:,3)=I_tot(:,:,3)*(1-DoP/DoP_infi)./Trans(:,:,3);

L_obj(:,:,1)=(L_obj(:,:,1)-min(min(L_obj(:,:,1))))/(max(max(L_obj(:,:,1)))-min(min(L_obj(:,:,1))));
L_obj(:,:,2)=(L_obj(:,:,2)-min(min(L_obj(:,:,2))))/(max(max(L_obj(:,:,2)))-min(min(L_obj(:,:,2))));
L_obj(:,:,3)=(L_obj(:,:,3)-min(min(L_obj(:,:,3))))/(max(max(L_obj(:,:,3)))-min(min(L_obj(:,:,3))));

A=(I2-I1)./DoP_infi;

figure(2)
imshow(L_obj);title('polar');

%% %%%%%%%%%%%%% Refine transimission %%%%%%%%
for i = 1:1:3
    refine_t(:,:,i) = RefineTrans(L_obj, Trans(:,:,i), 2, 0.5);
    
end

%% %%%%%%%%%%%% Refine radiance %%%%%%%%
 for i = 1:1:3
     L_opt(:,:,i) = max(max(L_obj(:,:,i)))*bilateralFilter(L_obj(:,:,i)/max(max(L_obj(:,:,i))),(1-refine_t(:,:,i)/max(max(refine_t(:,:,i)))));
 end


%% %%%%%%%%%%%% Final dehazing %%%%%%%%%%

figure(3)
imshow(L_opt);title('polar_opt');
imwrite(L_opt,'L01.jpg');




    
   




