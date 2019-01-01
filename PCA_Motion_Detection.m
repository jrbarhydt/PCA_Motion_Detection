% PCA Motion Detection
% Johnathon R Barhydt
clear all, close all, clc

%initialization variables
mode='';
a = 1/4;
N=2; %video set
lambda=0.005;

%pick and name video set
c1=strcat('cam1_',num2str(N),'.mat');
c2=strcat('cam2_',num2str(N),'.mat');
c3=strcat('cam3_',num2str(N),'.mat');

% get camera data
cam1 = importdata(c1);
cam2 = importdata(c2);
cam3 = importdata(c3);

%num_frames = min([size(cam1,4),size(cam2,4),size(cam3,4)]);
num_frames1= size(cam1,4);
num_frames2= size(cam2,4);
num_frames3= size(cam3,4);

%image size
im_h=uint32( size( cam1,1)*a);
im_w=uint32( size( cam1,2)*a);

%initialize frames
X1_images = zeros(im_h,im_w, num_frames1);
X2_images = zeros(im_h,im_w, num_frames2);
X3_images = zeros(im_h,im_w, num_frames3);
%populate frames
for i=1:num_frames1
   X1_images(:,:,i) = im2double( imresize(rgb2gray( cam1(:,:,:,i)),[im_h im_w]));
end
for i=1:num_frames2
   X2_images(:,:,i) = im2double( imresize(rgb2gray( cam2(:,:,:,i)),[im_h im_w]));
end
for i=1:num_frames3
   X3_images(:,:,i) = im2double( imresize(rgb2gray( cam3(:,:,:,i)),[im_h im_w]));
end

%% Run each video through vid2pos to collect x, y coordinates
[x1,y1]=vid2pos(X1_images,'noisy',20,lambda);
[x2,y2]=vid2pos(X2_images,'noisy',20,lambda);
[x3,y3]=vid2pos(X3_images,'noisy',20,lambda);
%% Crop and shift coords to same timebase
x1(1:20)=[];
x2(1:20)=[];
x3(1:20)=[];
y1(1:20)=[];
y2(1:20)=[];
y3(1:20)=[];

if var(x1)>var(y1)
    [~,ind]=max(x1(1:30));
else
    [~,ind]=max(y1(1:30));
end
x1(1:ind)=[];
y1(1:ind)=[];

if var(x2)>var(y2)
    [~,ind]=max(x2(1:30));
else
    [~,ind]=max(y2(1:30));
end
x2(1:ind)=[];
y2(1:ind)=[];

if var(x3)>var(y3)
    [~,ind]=max(x3(1:30));
else
    [~,ind]=max(y3(1:30));
end
x3(1:ind)=[];
y3(1:ind)=[];

trunc = min([size(x1,1),size(x2,1),size(x3,1)]);
x1(trunc:end)=[];
x2(trunc:end)=[];
x3(trunc:end)=[];
y1(trunc:end)=[];
y2(trunc:end)=[];
y3(trunc:end)=[];
%%
X=[x1,y2,x2,y2,x3,y3];
[L,S]=inexact_alm_rpca(X);
X=L;
[m,n]=size(X);
nuclear_norm=mean(X,2);
X=X-repmat(nuclear_norm,1,n); 
[u,s,v]=svd(X'/sqrt(n-1),'econ');
lambda=diag(s).^2;
l_norm = lambda/sum(lambda);

figure(1)
subplot(2,1,1)
bar(l_norm), title('Principle Component Strength - Trial 2'), xlabel('mode'),ylabel('strength in %')
subplot(2,1,2)
plot(v(:,1)), hold on, plot(v(:,2)), hold on, plot(v(:,3)), xlabel('time (frame#)'),ylabel('position (rel)')
legend('primary mode','secondary mode','tertiary mode')
%%
figure(1)
title('Quantitative Comparison')
subplot(6,1,1)
plot(x1), axis off
subplot(6,1,2)
plot(y1), axis off
subplot(6,1,3)
plot(x2), axis off
subplot(6,1,4)
plot(y2), axis off
subplot(6,1,5)
plot(x3), axis off
subplot(6,1,6)
plot(y3), axis off
