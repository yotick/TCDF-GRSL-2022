%% get detials by reducing scale
function [D, k] = get_D_red(Pc,img_mul)
global sensor ratio sate;

r=2; 
esp=10^-6; eps = 10^-6;
M =imresize(img_mul,size(Pc),'bicubic');   % 双三次插值算法    
[I,alpha] = get_I(sate,M, Pc); 
 
%%  degraded scale initial
pan_d = imresize(Pc,1/ratio);

ms_LP = MTF(img_mul,sensor,[],ratio);
ms_LP_d = imresize(ms_LP,1/ratio);
ms_up = imresize(ms_LP_d, size(pan_d));
Dm = img_mul-ms_up;
% pan_LP = MTF_PAN(Pc,sensor,ratio);
%     pan_d = pan_LP(3:ratio:end,3:ratio:end);


% ms_LP = MTF(img_mul,sensor,[],ratio);
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
end