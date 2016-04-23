clear all;
close all;
clc;

Img=imread('three.bmp');
 Img = rgb2gray(Img);
%  Img = Img(:,size(Img,1)/2-5:end);
%   ff = fspecial('laplacian');
%   Img = imfilter(Img,ff);

% % Img = rgb2gray(Img);
% figure(1); imshow(uint8(Img));
%      Img=imresize(Img,0.5);
% [m,n]=size(Img);
% %%- 将初始曲线C设置为圆
% 
%  Img = Img(10:end-25,10:end-10);
 c0=2;
 initialLSF=c0*ones(size(Img));
 % generate the initial region R0 as a rectangle
 initialLSF(10:size(Img,1)-10,3 :size(Img,2)-3)=-c0;  
 phi=initialLSF;
 

[nrow,ncol] =size(Img);
Img2=ones(nrow,ncol);
% get the size
[nrow,ncol] =size(Img);
% 
  ic=nrow/2;
  jc=ncol/2;
  r=ic/6;
 %  
   %phi_0 = sdf2circle(nrow,ncol,ic,jc,r);
   [X,Y] = meshgrid(1:ncol, 1:nrow);
  phi = sqrt((X-jc).^2+(Y-ic).^2)-r;  %初始化phi为SDF

delta_t = 0.1;
lambda_1=1;
lambda_2=1;
nu=0;
h = 1;
epsilon=1;
mu = 0.01*255*255; 

figure;%显示初始化图像 figure1
imshow(Img);colormap
hold on;

%plotLevelSet(phi_0,0,'g');
[c,h] = contour(phi,[0 0],'r'); 
hold off;
%  apparence(Img,phi,1.5,sqrt(50))
 for iter = 1:200
      figure(3)
        if mod(iter,50)==0
        imshow(Img);colormap
        hold on;
        [c,h] = contour(phi,[0 0],'r'); 
        hold off;
        pause(.5);
        end
      phi = apparent(Img,phi,1,sqrt(50)); 
    
%         figure(3)
%         if mod(iter,1)==0
%         imshow(Img);colormap
%         hold on;
%         [c,h] = contour(phi,[0 0],'r'); 
%         hold off;
%         pause(5)
%         end
%      
 end
 figure(4);
 phi_inverse = phi(end:-1:1,end:-1:1);
 [c,h] = contour(phi_inverse,[0,0],'r');


  
        


