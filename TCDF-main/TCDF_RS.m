clc; clear all; close all;
%% only for plaidas
addpath(genpath('..\shared'))

global  ratio file_path_rgb256;
global count time curr_d gt256;
global sensor sate num mu L Qblocks_size flag_cut_bounds dim_cut thvalues;
% global eps;

curr_d = pwd;
sate = 'ik';   % geo,ik,pl, qb, wv2，wv3

method = 'TCGF';

%% start initialization
initialize_sate_RS_MTF();

%% start process
% for  num= 1 : length(file_path_rgb256)
for  num= 10:10
    count = 1 + count
    %% read images and preprocess, the users should defined by themselfs.
    [rgb256, gt256, mul64, pan256 ] = read_image_sate_RS(num);
    
    %% tic start
    tic;
    ORM = gt256;
    img_mul = im2double(mul64);
    P = im2double(pan256);
    %     P_L = MTF_PAN(P,sensor,ratio);
    Pd = imresize(P,1/4);
    
    %% 加权平均求I分量
    [m,n] = size(P);
    M =imresize(img_mul,size(P),'bicubic');   % 双三次插值算法
    img_mul_L = MTF(img_mul,sensor,0,ratio);
    [I,alpha] = get_I(sate,M, P);

    %% ************ function of getting Pc*******

    [Pc, De, H] = get_Pnew(P,I);
    [D_new] = get_D_red2(Pc,img_mul); 

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
            g = 1.21;
        case 'qb'
            g = 0.80;
        case 'geo'
            g = 1.35;
    end  
     F_final2 = zeros(size(M)); 
    for i=1:4

        F_final2(:,:,i) = M(:,:,i)+ g*M_rate(:,:,i).* (D_new(:,:,i));
    end
    

    %% end and get RGB
    toc
    time(num) = toc;    

    for i = 1:3
        F_rgb2(:,:,i) =  F_final2(:,:,i);
    end

    %% save images
    save_sate_RS(F_final2, M, F_rgb2, num);
end

%% show image
% saveData_RS(sig_all,co_PLI, co_PI);
Eval = Evaluation4(ORM,pan256,uint8(F_final2*255))
[Q_avg_Segm, SAM_Segm, ERGAS_Segm, SCC_GT_Segm, Q_Segm] = indexes_evaluation(...
    uint8(F_final2*255),ORM,ratio,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues);
Eval2 = [Q_avg_Segm, SAM_Segm, ERGAS_Segm, SCC_GT_Segm, Q_Segm]

% figure,
% imshow(F_rgb);
% title('After NVDI');
figure,
imshow(F_rgb2);
title('Before NVDI');
figure,
imshow(rgb256);
title('GT');
% figure, imshow(rgb256)
cd(curr_d);
%% write result
T = quality_eval_sate_RS();
cd(curr_d);
writetable(T,strcat('final_result/quality_',method,'_4C_',sate,'_RS_new4.csv'),'WriteRowNames',true);

