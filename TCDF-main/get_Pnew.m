function [IH2, De, H] = get_Pnew(P,I)
global  sate;

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