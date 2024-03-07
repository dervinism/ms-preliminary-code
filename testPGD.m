ContinuesPhimatrix = [0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 0 0 0 0 0 0;
                      0 0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 0 0 0 0 0;
                      0 0 0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 0 0 0 0;
                      0 0 0 0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 0 0 0;
                      0 0 0 0 0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 0 0;
                      0 0 0 0 0 0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 0;
                      0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5];

[Gx,Gy] = gradient(flip(ContinuesPhimatrix));

Mag_of_MeanGrad = sqrt(sum(nansum(Gx))^2+sum(nansum(Gy))^2);
Mean_of_MagGrad = sum(nansum(sqrt(Gx.^2+Gy.^2)));
%PGD = Mag_of_MeanGrad/Mean_of_MagGrad
(sum(nansum(Gx.*Gy)))/sqrt(sum(nansum(Gx.^2))*sum(nansum(Gy.^2)))



ContinuesPhimatrix = [0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 0 0 0 0 0 0 0 0 0 0 0 0;
                      0 0 0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 0 0 0 0 0 0 0 0 0 0;
                      0 0 0 0 0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 0 0 0 0 0 0 0 0;
                      0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 0 0 0 0 0 0;
                      0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 0 0 0 0;
                      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 0 0;
                      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5];

[Gx,Gy] = gradient(flip(ContinuesPhimatrix));

Mag_of_MeanGrad = sqrt(sum(nansum(Gx))^2+sum(nansum(Gy))^2);
Mean_of_MagGrad = sum(nansum(sqrt(Gx.^2+Gy.^2)));
%PGD = Mag_of_MeanGrad/Mean_of_MagGrad
(sum(nansum(Gx.*Gy)))/sqrt(sum(nansum(Gx.^2))*sum(nansum(Gy.^2)))



ContinuesPhimatrix = [0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 1 0 0 0 0 0 0 0 0 0 0 0;
                      0 0 0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 1 0 0 0 0 0 0 0 0 0;
                      0 0 0 0 0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 1 0 0 0 0 0 0 0;
                      0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 1 0 0 0 0 0;
                      0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 1 0 0 0;
                      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 1 0;
                      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5];

[Gx,Gy] = gradient(flip(ContinuesPhimatrix));

Mag_of_MeanGrad = sqrt(sum(nansum(Gx))^2+sum(nansum(Gy))^2);
Mean_of_MagGrad = sum(nansum(sqrt(Gx.^2+Gy.^2)));
%PGD = Mag_of_MeanGrad/Mean_of_MagGrad
(sum(nansum(Gx.*Gy)))/sqrt(sum(nansum(Gx.^2))*sum(nansum(Gy.^2)))



ContinuesPhimatrix = flipud([0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 0 0 0 0 0 0 0 0 0 0 0 0;
                      0 0 0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 0 0 0 0 0 0 0 0 0 0;
                      0 0 0 0 0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 0 0 0 0 0 0 0 0;
                      0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 0 0 0 0 0 0;
                      0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 0 0 0 0;
                      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 0 0;
                      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5;
                      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5]);

[Gx,Gy] = gradient(flip(ContinuesPhimatrix));

Mag_of_MeanGrad = sqrt(sum(nansum(Gx))^2+sum(nansum(Gy))^2);
Mean_of_MagGrad = sum(nansum(sqrt(Gx.^2+Gy.^2)));
%PGD = Mag_of_MeanGrad/Mean_of_MagGrad
(sum(nansum(Gx.*Gy)))/sqrt(sum(nansum(Gx.^2))*sum(nansum(Gy.^2)))



% ContinuesPhimatrix = [0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 0 0 0 0 0 0 0 0 0 0 0 0;
%                       0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 0 0 0 0 0 0 0 0 0 0 0 0;
%                       0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 0 0 0 0 0 0 0 0 0 0 0 0;
%                       0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 0 0 0 0 0 0 0 0 0 0 0 0;
%                       0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 0 0 0 0 0 0 0 0 0 0 0 0;
%                       0 0 0 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 0 0 0 0 0 0 0 0 0 0 0 0;
%                       0 0 0 0 0 0 0 1 1 2 2 2 3 3 3 4 4 4 5 5 5 0 0 0 0 0 0 0 0 0 0 0 0];
% 
% [Gx,Gy] = gradient(flip(ContinuesPhimatrix));
% 
% Mag_of_MeanGrad = sqrt(sum(nansum(Gx))^2+sum(nansum(Gy))^2);
% Mean_of_MagGrad = sum(nansum(sqrt(Gx.^2+Gy.^2)));
% %PGD = Mag_of_MeanGrad/Mean_of_MagGrad
% (sum(nansum(Gx.*Gy)))/sqrt(sum(nansum(Gx.^2))*sum(nansum(Gy.^2)))



K = NaN;
cntr = -1;
while isnan(K)
  cntr = cntr+1;
  %K = (ContinuesPhimatrix(floor(end/2)+cntr,floor(end/2)+cntr)/(2*pi));
  K = (ContinuesPhimatrix(min([end ceil(end/2)+cntr]),min([end ceil(end/2)+cntr]))/(2*pi));
end

if sum(nansum(Gx(:,:)))<0
  NetGradDir =  pi+atan(sum(nansum(Gy(:,:)))/sum(nansum(Gx(:,:))))
else
  NetGradDir=  atan(sum(nansum(Gy(:,:)))/sum(nansum(Gx(:,:))))
end