function [x,y] = vid2pos(A,mode,slices,lambda,fail)
%This program takes in a matrix of video frames and attempts to find
%the location of a single moving object as it moves from frame to  frame
%
%There are four 'modes' to compare performance of different techniques
%including:         mode                      string
%                 -----------------------------------
%                 Naive BG subtraction        'naive'
%                 SVD low-rank subtraction    'clean'
%                 Robust ALM PCA Sparse BG    'noisy'
%                 PCA Image Projection        DEFAULT
%
% Additionally, there are two parameters, slices and lambda
% slices: accepts an integer number to slice the video into segments
% lambda: accepts a value<1 indicating sparsity of rPCA method
%            see the following website/file: 
%http://perception.csl.illinois.edu/matrix-rank/Files/inexact_alm_rpca.zip
%   fail: back-up mode if windowing is ineffective, set to 'yes' or blank

    %default arguments
    if nargin < 2
       mode = 'clean'; 
    end
    if nargin < 3
       slices = 2; 
    end
    if nargin < 4
        lambda = 0.05;
    end
    if nargin < 5
        fail='';
    end
    
    % initialization variables
    num_frames=size(A,3);
    im_h=size(A,1);
    im_w=size(A,2);
    clip_length = uint32( num_frames/slices);
    x=zeros(num_frames,1);
    y=zeros(num_frames,1);
    
    %construct images into stacked matrix, each column is one frame
    for i=1:num_frames
        X(:,i) = reshape( A(:,:,i), im_w*im_h, 1);
    end
    
    %Below are the operation 'modes'
    
    %Robust PCA via Augmented Lagrange Multiplier
    if isequal(mode,'noisy')
        [~,S]=inexact_alm_rpca(X,lambda);
    
    %SVD Low-Rank Approxmiation
    elseif isequal(mode,'clean')
        [u,s,v]=svd(X,'econ');
        q=10;
        L=u(:,1:q)*s(1:q,1:q)*v(:,1:q).';
        S=X-L;

    %Simple Background Frame Subtraction
    elseif isequal(mode,'naive')
        S=X(:,1:end-1)-X(:,2:end);
        num_frames=num_frames-1; 

    %Standard PCA method
    else
        [m,n]=size(X);
        nuclear_norm=mean(X,2);
        Xr=X-repmat(nuclear_norm,1,n); 
        [u,s,v]=svd(Xr'/sqrt(n-1),'econ');
        Y=Xr*u';L=u(:,1:3)*s(1:3,1:3)*v(:,1:3).';
        for i=1:num_frames
            [x(i+1),y(i+1)]=locator(reshape(Y(:,i),im_h,im_w),'point',[x(i), y(i)]);
        end
        return
    end

    %Run in simple,clean mode if windowing fails
    if isequal(fail,'yes')
        for i=1:num_frames
            frame=reshape(S(:,i),im_h,im_w);
            frame=frame.*frame;
            [x(i),y(i)]=locator(frame.*(frame>0));
        end
        return
    end
    
    % determine iterations for various slice sizes
    for j=1:slices
        if j==slices&&slices>1
            if mod(num_frames,clip_length)==0
                iters=clip_length;
            else
                iters=mod(num_frames,clip_length);
            end
        else
            iters=clip_length;
        end

        %Initialize training values for window location
        x_train = [];
        y_train = [];
        for i=1:iters
            k=(j-1)*clip_length+i;
            [x_train(i),y_train(i)]=locator(reshape(S(:,k),im_h,im_w));
        end

        %Remove null entries and take training set mean location
        x_train(x_train==0) = [];
        y_train(y_train==0) = [];
        x(k)=mean(x_train);y(k)= mean(y_train);

        %Run locator in 'point' mode, performing moving-window capability
        for i=1:iters
            k=(j-1)*clip_length+i;
            [x(k+1),y(k+1)]=locator(reshape(S(:,k),im_h,im_w),'point',[x(k), y(k)]);
        end       
    end

