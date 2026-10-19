%% This program is used for downsampling PSF
% last update:0603. Auth: MW
clc;clear;
close all;
Nnum = 15;

% psfpath = strcat('D:\RUSH3Dproject\RUSH3Dpsfblock\20220818\ROI',num2str(id),'_5.761x_lambda_525_axial_-400um_to_400um_zstep_8_NA_0.5_fml_393.3um_OSR_5');

psfpath = 'D:\RUSH3Dproject\RUSH3Dpsfsim\20240329\r1_5.761x_lambda_525_axial_-300um_to_300um_zstep_6_NA_0.5_fml_393.3um_ftl_265mm_OSR_5';
psfmatpath = strcat(psfpath,'\psf_zmat');
fileall = dir(fullfile(psfmatpath,'*.mat'));
files_name = sort_nat({fileall.name});
downsampling_rate = 3 ;
savefilepath = strcat(psfpath,'_ds_',num2str(downsampling_rate),'\psf_zmat');
savefilepath_pic = strcat(psfpath,'_ds_',num2str(downsampling_rate),'\psf_z');
mkdir(savefilepath);
mkdir(savefilepath_pic);
psf_z = load(strcat(psfmatpath,'\\',files_name{1}));
if isstruct(psf_z )
    psf_z  = psf_z.psf_z;
end
size_new = floor(size(psf_z,1)/2/downsampling_rate)*2+1;
for z = 1 : length(files_name)
    psf_z = load(strcat(psfmatpath,'\\',files_name{z}));
    % handle the structure
    if isstruct(psf_z )
        psf_z  = psf_z.psf_z;
    end
    psf_z = imresize(psf_z,[size_new,size_new]);
    
    Output=uint16(gather(psf_z)./max(psf_z(:)).*65535);
    for idu=1:Nnum
        for idv = 1:Nnum
            ii = (idu-1)*Nnum+idv;
            if ii==1
                imwrite(uint16(squeeze(Output(:,:,idu,idv))),strcat(savefilepath_pic,'\psf_sim',num2str(z),'.tif'));
            else
                imwrite(uint16(squeeze(Output(:,:,idu,idv))),strcat(savefilepath_pic,'\psf_sim',num2str(z),'.tif'),'WriteMode', 'append');
            end
        end
    end
    
    save(strcat(savefilepath,'\psf_sim',num2str(z),'.mat'),'psf_z','-v7.3');
    disp([num2str(z),' has been ds...']);
end

