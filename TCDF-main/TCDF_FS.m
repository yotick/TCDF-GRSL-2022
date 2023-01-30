clc; clear all; close all;
%% only for plaidas 
addpath(genpath('E:\remote sense image fusion\shared'))

global   thvalues  ratio L; 
global   im_tag sensor;
global   file_path_rgb_noR; 
global  count time num ;
global  sate curr_d;

curr_d = pwd;
sate = 'ik';   % geo,ik,pl, qb, wv2，wv3

method = 'TCGF';
%% start initialization
initialize_sate_FS();

%% start process
% for  num= 1 : length(file_path_rgb_noR)t
 for  num= 1:5
     count = 1 + count
    %% read images and preprocess
    [mul_noR, pan_noR ] = read_image_sate_FS(num);
    
    %% tic start
    tic()   
%     ORM = gt256;
    img_mul = im2double(mul_noR);
    P = im2double(pan_noR);    
    
       %% 加权平均求I分量
    [m,n] = size(P);
    M =imresize(img_mul,size(P),'bicubic');   % 双三次插值算法    
 
    img_mul_L = MTF(img_mul,sensor,0,ratio);
    [I,alpha] = get_I(sate,M, P);

    %% ************ function of getting Pc*******

    [Pc, De, H] = get_Pnew(P,I);
    [D] = get_D_red2(Pc,img_mul); 
    D_new = D;

%% *********** new code

    
    %% original
    M_sum = M(:,:,1)+M(:,:,2)+M(:,:,3)+M(:,:,4);
    M_sum(M_sum==0)=eps;
    M_rate = zeros(size(M));
    for i=1:4
        M_rate(:,:,i) = 4*M(:,:,i)./M_sum;
    end
    

    %% switch
    
    switch sate
        case 'ik'
            g = 1;
        case 'pl'
            g = 1;
        case 'wv2'
            g = 1.21;
        case 'wv3'
            g = 1;
        case 'qb'
            g = 0.80;
        case 'geo'
            g = 1.35;
    end    
    

     F_final = zeros(size(M)); 
    for i=1:4

        F_final(:,:,i) = M(:,:,i)+ g*M_rate(:,:,i).* (D_new(:,:,i));
    end
    

    %% end and get RGB
    toc
    time(num) = toc;
    

    for i = 1:3
        F_rgb(:,:,i) =  F_final(:,:,i);
    end
     %% save images
        save_sate_FS(F_final, M, F_rgb, num);
end

 %% show image
figure, imshow(F_rgb);
[Dl,Ds,QNR_index,SAM_index,sCC] = indexes_evaluation_FS(F_final,img_mul,P,...
    L,thvalues,M,sensor,im_tag,ratio);
% [D_lambda,D_S,QNR_index,SAM_index,sCC] = indexes_evaluation_FS(Fused,I_MS_LR,I_PAN);
Eval = [Dl,Ds,QNR_index,SAM_index,sCC]

%% write result
 T = quality_eval_sate_FS();
cd(curr_d);
writetable(T,strcat('final_result/quality_',method,'_4C_',sate,'_FS2.csv'),'WriteRowNames',true);
