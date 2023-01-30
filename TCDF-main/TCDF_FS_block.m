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
 for  num= 13:13
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

    
        %% block processing
    S_block = 256;
    
    in3 = cat(3,P,I);
    fun_inj = @(bs)get_Pnew(bs.data);
    Pc = blockproc(in3,[S_block, S_block],fun_inj);
    img_up = imresize(img_mul,4, 'nearest');
    in3 = cat(3,Pc,img_up);
    fun_inj = @(bs)get_D_red2(bs.data);
    D = blockproc(in3,[S_block, S_block],fun_inj);
 
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
writetable(T,strcat('final_result/quality_',method,'_4C_',sate,'_FS_B.csv'),'WriteRowNames',true);

%% ************* function ***************************************
function [IH2] = get_Pnew(in3)
global  sate;

P = in3(:,:,1);
I= in3(:,:,2);
% P=(P-mean(P(:)))*std2(I)/std(P(:)) + mean2(I);   % histogram matching
[H, sig, cor,PL] = get_gau_H_RS(P,I);
switch sate
    case 'ik'
        u = 85;
    case 'pl'
        u = 50;
    case 'wv2'
        %  u = get_u(P,H,I);
        u = 50;
    case 'wv3'
        u = 50;
        %            u = get_u(P,H,I);
    case 'qb'
        u = 5;
    case 'geo'
        u = 50;
end
%     P_L = ifft2(fft2(P).*H);
%     imshow(P_L)
%     figure, imshow(I);

lap = fspecial('laplacian',0);
lap_f = freqz2(lap,size(P)); % laplase matrix
lap_f = ifftshift(lap_f);
%% original
A = H.*fft2(I) + u*lap_f.*lap_f.* fft2(P);
B = H.*H + u*lap_f.*lap_f;

C = A./B;
IH2 = real(ifft2(C));
De = IH2 - I;
end

%% *********  function2 *************************************

function [Df, k] = get_D_red2(in3)
global sensor ratio sate;

Pc = in3(:,:,1);
img_up = in3(:,:,2:5);
img_mul = imresize(img_up,1/4,'nearest');
r=2;
esp=10^-6;
M =imresize(img_mul,size(Pc),'bicubic');   % 双三次插值算法
[I,alpha] = get_I(sate,M, Pc);

ML =  MTF(M,sensor,[],ratio);
D_M = M - ML;
%%  degraded scale initial
pan_d = imresize(Pc,1/ratio);

ms_LP = MTF(img_mul,sensor,[],ratio);
ms_LP_d = imresize(ms_LP,1/ratio);
ms_up = imresize(ms_LP_d, size(pan_d));
Dm = img_mul-ms_up;
% pan_LP = MTF_PAN(Pc,sensor,ratio);
%     pan_d = pan_LP(3:ratio:end,3:ratio:end);


ms_LL = MTF(ms_LP,sensor,[],ratio);
Dmm = ms_LP-ms_LL;
Dmm_ave = 1/4*(Dmm(:,:,1)+Dmm(:,:,2)+Dmm(:,:,3)+Dmm(:,:,4));

% ms_LP_d = imresize(ms_LP,1/ratio );
% ms_up = imresize(ms_LP_d, ratio);

%      pan_d = histog( ms_up, pan_d);

%     Dm = img_mul-ms_up;
I_orig = alpha(1)*img_mul(:,:,1) + alpha(2)*img_mul(:,:,2)+alpha(3)*img_mul(:,:,3)+alpha(4)*img_mul(:,:,4);
%     ms_ex = ms_up;
I_d = alpha(1)*ms_up(:,:,1) + alpha(2)*ms_up(:,:,2)+alpha(3)*ms_up(:,:,3)+alpha(4)*ms_up(:,:,4);

[H_d,sig_d,cor_old_d, P_L_d] = get_gau_H_RS(pan_d, I_d);
p_L1 = ifft2(fft2(pan_d).*H_d);
% p_L1 = imguidedfilter(pan_d, I_d);
% p_L1 = guidedfilter(pan_d,I_d,r,esp);
% p_L1 = gradient_guidedfilter(pan_d,I_d, eps);
% H2 = get_H_MTF_P(pan_d,sensor,ratio);    %% 转换为FFT频率域可行 in IKONO is OK,
% p_L2= ifft2(fft2(pan_d).*H2);

ms_ex(:,:,1) = I_d;
ms_ex(:,:,2)= p_L1;
% ms_ex(:,:,3)= p_L2;
%     ms_ex(:,:,4)= p_L3;

%% full scale  initial

[H,sig,cor_old, P_L_ori] = get_gau_H_RS(Pc, I);
P_L1 = ifft2(fft2(Pc).*H);
% P_L1 = gradient_guidedfilter(Pc, I, eps);
% P_L1 = imguidedfilter(Pc, I);
% P_L1=guidedfilter(Pc,I,r,esp);
% H4 = get_H_MTF_P(Pc,sensor,ratio);    %% 转换为FFT频率域可行 in IKONO is OK,
% P_L2 = ifft2(fft2(Pc).*H4);

%% start for circle
for i = 1:4
    %% degrade scale
    
    beta(i,:) = impGradDes(ms_ex, pan_d - Dm(:,:,i));  % 4*5
    %       beta(i,:) = estimation_alpha(ms_ex,pan_d - Dm(:,:,i),'global');
    %     p_LP_d2(:,:,i) = beta(i,1)*I_d + beta(i,2)*p_L1 + beta(i,3)*p_L2 ;
    p_LP_d2(:,:,i) = beta(i,1)*I_d + beta(i,2)*p_L1 ;
    Dd(:,:,i) =pan_d -p_LP_d2(:,:,i);
    %       p_LP_d2(:,:,i) = beta(i,1)*I_d + beta(i,2)*p_L2 ;
    %% end
    
    
    %     k(i) = impGradDes( Dd(:,:,i), Dm(:,:,i));
    k(i) = impGradDes(Dm(:,:,i), Dd(:,:,i));
    %     k2(i) = mean2( abs(Dm(:,:,i)))/mean2(abs(Dd(:,:,i)));
    %         Dp_d(:,:,i) =  pan_d - p_LP_d2(:,:,i);
    
    %% full scale
    
    %     P_LP(:,:,i) = beta(i,1)*I+beta(i,2)*P_L1+beta(i,3)*P_L2 ;
    P_LP(:,:,i) = beta(i,1)*I+beta(i,2)*P_L1;
    %         P_LP(:,:,i) = beta(i,1)*I+beta(i,2)*P_L2;
    %% end
    D(:,:,i) = Pc-P_LP(:,:,i);
%     D(:,:,i) = Pc-I;  % ablation
end

Dd_ave = 1/4*(Dd(:,:,1)+Dd(:,:,2)+Dd(:,:,3)+Dd(:,:,4));
% t = mean2(abs(Dd_ave))/mean2(abs(Dmm_ave));
% t2 = impGradDes(Dmm_ave, Dd_ave);
% Dmm = t2*Dmm;
Df = zeros(size(M));
tao = zeros(4,2);
for i = 1:4
%     t2(i) = impGradDes(Dmm(:,:,i), Dd(:,:,i));
%     Dmm(:,:,i) = t2(i)*Dmm(:,:,i);
    D1_ex(:,:,1) = Dd(:,:,i);
    D1_ex(:,:,2) = Dmm(:,:,i);
    tao(i,:) = impGradDes(D1_ex, Dm(:,:,i));
    if strcmp(sate,'wv3')==1 || strcmp(sate,'wv2')==1|| strcmp(sate,'geo')==1  
         tao(i,1) = max(tao(i,1), 0.85);
    end   
    Df(:,:,i) = tao(i,1)*D(:,:,i) + tao(i,2)*D_M(:,:,i);
end
% Df = D;  % ablation
end
    
