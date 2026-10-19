clc;clear;
close all;

iter = 1;
Nnum = 15;
psfpath = 'D:\RUSH3Dproject\RUSH3Dpsf\46.7x_lambda525_axial_-18um_to_18um_zstep_0.4_NA_1.1_fml_2100um_ftl_180mm_OSR_5\psf_zmat';
psfabpath = 'D:\RUSH3Dproject\RUSH3Dpsf\SA-2_46.7x_lambda_525_axial_-18um_to_18um_zstep_0.4_NA_1.1_fml_2100um_ftl_180mm_OSR_5\psf_zmat';
load(sprintf('%s\\psf_sim60.mat',psfpath)); %psf_z
PSFsimu = psf_z;
load(sprintf('%s\\psf_sim60.mat',psfabpath));
PSFaber = psf_z;


map_wavshape = zeros(Nnum,Nnum,2);
[Sx,Sy] = meshgrid(-fix(Nnum/2):fix(Nnum/2),-fix(Nnum/2):fix(Nnum/2));
mask = (Sx.^2+Sy.^2) <= 25;
for u = 1:Nnum
    for v = 1:Nnum
        PSFsimu_uv = PSFsimu(:,:,u,v);
        PSFreal_uv = PSFaber(:,:,u,v);
        corr_map=normxcorr2(PSFsimu_uv,PSFreal_uv);
        [shift_a,shift_b]=find(corr_map==max(corr_map(:)));
        
        map_wavshape(u,v,1)=shift_a(1)-size(PSFsimu_uv,1);
        map_wavshape(u,v,2)=shift_b(1)-size(PSFsimu_uv,2);
    end
    disp(u);
end
%%
map_wavshape(:,:,1) = (map_wavshape(:,:,1) - map_wavshape(ceil(Nnum/2),ceil(Nnum/2),1)).*mask;
map_wavshape(:,:,2) = (map_wavshape(:,:,2) - map_wavshape(ceil(Nnum/2),ceil(Nnum/2),2)).*mask;
%%
savepath = 'D:\03project\01RUSH3D\RUSH3Dproject\RUSH3Dresult\psfaberation_gt\map';
mkdir(sprintf('%s\\iter_%d',savepath,iter));
save(sprintf('%s\\iter_%d\\map_wavshape%d.mat',savepath,iter,iter),'map_wavshape');