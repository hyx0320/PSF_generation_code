function aberration = create_aberration(D,ratio)

x = linspace(-1,1,D);%%采样率

[kx, ky] = meshgrid(x,x); %%离散化点格
[theta,r] = cart2pol(kx,ky);
idx = r <= 1;
Dic = [0,0;1,-1;1,1;2,-2;2,0;2,2;3,-3; 3,-1 ;3,1; 3,3; 4,-4;4,-2;4,0;4,2;4,4;5,-5; 5,-3; 5,-1 ;5,1; 5,3; 5,5; 6,-6 ;6,-4; 6,-2; 6,0; 6,2; 6,4; 6,6;...
    7,-7; 7,-5 ; 7,-3 ; 7,-1 ; 7,1 ; 7,3 ; 7,5 ; 7,7 ; 8,-8 ; 8,-6 ; 8,-4 ; 8,-2 ; 8,0; 8,2; 8,4; 8,6; 8,8]; %%阶数
aberration = zeros(size(kx));

% for i = 7:1:45
%     weight = rand(1);
%     if i == 7
%         aberration(idx) = weight*zernfun(Dic(i,1),Dic(i,2),r(idx),theta(idx));
%     else
%         aberration(idx) =aberration(idx) + weight*zernfun(Dic(i,1),Dic(i,2),r(idx),theta(idx));
%     end
% end

aberration(idx) = zernfun(Dic(13,1),Dic(13,2),r(idx),theta(idx));%+2*zernfun(Dic(13,1),Dic(13,2),r(idx),theta(idx));

rms=calcRMS(aberration);
aberration=aberration./rms.*2*pi.*ratio;% 一个2*pi，对应一个lambda

aberration(isnan(aberration)) = 0;

end