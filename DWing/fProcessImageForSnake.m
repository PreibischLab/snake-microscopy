

function I = fProcessImageForSnake(I, r, sigma, gHood)

G = fspecial('gaussian',gHood,sigma);
I = imfilter(I,G,'same');
I = imbothat(I, strel('disk', r));
I = f_LocalNormalization(I, 5, 20);
% I = I-3 *IB;

end

function I = f_LocalNormalization(I, sigma1, sigma2)
eps = 1e-1;
Size1 = 2*ceil(-norminv(eps/2,0, sigma1))+1;
Size2 = 2*ceil(-norminv(eps/2,0, sigma2))+1;
G1 = fspecial('gaussian', Size1,sigma1);
G2 = fspecial('gaussian', Size2,sigma2);
IdiffMean = I - imfilter(I,G1);
Istd = sqrt(imfilter(IdiffMean.^2,G2));
Istd(Istd==0) = 1;
I = IdiffMean./Istd;

end