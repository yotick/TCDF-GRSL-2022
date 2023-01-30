function a = get_co_grad(M_ori,Pd)

lap = fspecial('laplacian',0);
 grad_P = imfilter(Pd,lap);
%  imshow(grad_P);
 [m,n,z] = size(M_ori);
 grad_M = zeros(m,n,z);
 a = zeros(1, size(M_ori,3));
 for i = 1:size(M_ori,3)
     grad_M(:,:,i) = imfilter(M_ori(:,:,i),lap);
     a(i) = impGradDes(grad_M(:,:,i), grad_P);
%      a(i) = impGradDes(grad_P, grad_M(:,:,i) );
 end
 
end