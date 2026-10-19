clc;clear;
close all;
addpath('D:\Deepwonder2.0\pipeline_experiment1\Realign_timelaps');
view_config_file = 'view_circle5';
psfpath = 'D:\RUSH3Dproject\RUSH3Dpsf\4.561x_lambda525_axial_-300um_to_300um_zstep_6_NA_0.4_fml_393.3um_ftl_200mm_OSR_5';
z_pos = (-300:6:300)';
psfmatpath = strcat(psfpath,'\psf_zmat');
fileall = dir(fullfile(psfmatpath,'*.mat'));
files_name = sort_nat({fileall.name});
Nnum = 15;

savefilepath = strcat(psfpath,'_pt','\psf_zmat');
savefilepath_pic = strcat(psfpath,'_pt','\psf_z');
mkdir(savefilepath);
mkdir(savefilepath_pic);
load(strcat(psfmatpath,'\\',files_name{1}));
[h,w] = size(psf_z,[1,2]);

for z = 1 : length(files_name)
    load(strcat(psfmatpath,'\\',files_name{z}));
    psf_pt_pos = zeros(Nnum*Nnum,2);
    for idu=1:Nnum
        for idv = 1:Nnum
            point_mat = zeros(h,w);
            psf_zuv = psf_z(:,:,idu,idv);
            [max_v, ind] = max(psf_zuv(:));
            point_mat(ind) = 1;
            [x,y] = ind2sub([h,w],ind);
            point_mat = reshape(point_mat,[h,w]);
            ii = (idu-1)*Nnum+idv;
            psf_pt_pos(ii,:) = [x,y];
            if ii==1
                imwrite(im2uint16(point_mat),strcat(savefilepath_pic,'\psf_pt',num2str(z),'.tif'));
            else
                imwrite(im2uint16(point_mat),strcat(savefilepath_pic,'\psf_pt',num2str(z),'.tif'),'WriteMode', 'append');
            end
        end
    end
    
    save(strcat(savefilepath,'\psf_pt',num2str(z),'.mat'),'psf_pt_pos','-v7.3');
    disp([num2str(z),' processed...']);
end
%% load view configuration
View_config = imread(sprintf('%s.conf.sk.png',view_config_file));
v_list = zeros(sum(View_config(:)),2);
count = 0;
for i = 1 : Nnum*Nnum
    [u, v] = ind2sub([Nnum,Nnum],i);
    if View_config(v,u) > 0
        count = count + 1;
        v_list(count, :) = [v, u];
    end
end

%% calculate distance
count = 0;
psffit_matrix = zeros(size(v_list, 1)*(size(v_list, 1)-1)/2,6);
for i = 1 :size(v_list, 1) - 1
    for j = i + 1 : size(v_list, 1)
        count = count + 1;
        vindi = sub2ind([Nnum, Nnum], v_list(i,2),v_list(i,1));
        vindj = sub2ind([Nnum, Nnum], v_list(j,2),v_list(j,1));
        diffpos = zeros(length(files_name),2);
        for z = 1 : length(files_name)
            load(sprintf('%s\\psf_pt%d.mat',savefilepath,z),'psf_pt_pos');
            diffpos(z,:) = psf_pt_pos(vindi,:) - psf_pt_pos(vindj,:);
        end
        [p1,error1] = polyfit1_rmoutliers(z_pos,diffpos(:,1),2);
        [p2,error2] = polyfit1_rmoutliers(z_pos,diffpos(:,2),2);
        psffit_matrix(count,1) = vindi;
        psffit_matrix(count,2) = vindj;
        psffit_matrix(count,3:4) = p1;
        psffit_matrix(count,5:6) = p2;
        fprintf('View %d and View %d err: %.4f, %.4f...\n', vindi, vindj, error1, error2);
    end
end
save(sprintf('%s_pt\\psffit_matrix_%s.mat',psfpath, view_config_file),'psffit_matrix','-v7.3');