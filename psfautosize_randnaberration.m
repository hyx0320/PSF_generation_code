%% This program is used for calculate ideal PSF with no abberation
% Author: MW. Last update : 2021/05/04 
% psf: different view
% H: different position
%% Start
clc;clear;
close all;
gpuDevice(1);

%% parameter setting

% ------------- Basic parameters of sLFM ------------- %

M =         5.761; % M = NA*2*fml/MLPitch;% magnification || formula(15) ?
Nnum =      15; % 15 virtual pixels per microlens
n =         1; % refractive index
NA =        0.5; % NA of microlens and objective lens is 0.5. (NA of microlens and objective should be identical)
MLPitch =   69*1e-6; % The pitch size of microlens is 69 um
fml =       393.3*1e-6; % microlens focal length 193 um
ftl =       265e-3; % Focal length of tube lens is 165 mm 
OSR =       5; % calculate 5 times per pixel
lambda =    525*1e-9; %  light wave length is 525 nm

% -------------- Axial range and sampling rate ---------------- %

zmax_u = 300;
zmin_u = -300;
delta = 0.0001;
zmax = (zmax_u + delta)*1e-6;
zmin = (zmin_u + delta)*1e-6;
zspacing = 2e-6;%500e-6;

% --------------- Threshold ------------------ %

eqtol = 1e-10; % threshold no meaning = 0
tol = 0.0005; % threshold
Theshold  = 1e-3;

% ------------- Calculated parameters -------------- %

k = 2*pi*n/lambda; % k = 2pi/lambda x n (k = n x k0)
k0 = 2*pi*1/lambda; % k0 = 2pi/lambda
d = fml; % optical distance between the microlens and the sensor
fobj = ftl/M; % Focal length of objective lens 
fnum_obj = M/(2*NA); % f-number of objective lens (imaging-side) || not used in this code
fnum_ml = fml/MLPitch; % f-number of micrl lens || not used in this code
pixelPitch = MLPitch/Nnum; % pitch of virtual pixels *** important

if mod(Nnum,2)==0
    error('Nnum should be an odd number');
end

% --------------------- aberration ------------------------ %

aberration_flag = 1;

%r = [2,2.5,1.5,3]; % what meaning?
iter = 4;
%% Save PSF path
for ratio = 1.8
  
psf_folder1 = 'D:\03project\01RUSH3D\RUSH3Dproject\RUSH3Dpsf'; 
psf_folder_parameter = strcat('SA',num2str(ratio),'_',num2str(M),'x_lambda_',num2str(lambda*1e9),'_axial_',num2str(zmin_u),'um_to_',num2str(zmax_u), ...
                    'um_zstep_',num2str(zspacing*1e6),'_NA_',num2str(NA),'_fml_', ...
                    num2str(fml*1e6),'um_ftl_',num2str(ftl*1e3),'mm_OSR_',num2str(OSR));

psf_folder = strcat(psf_folder1, '\\', psf_folder_parameter);
psf_folder_H = strcat(psf_folder,'\\', 'H_z');
psf_folder_psf = strcat(psf_folder,'\\', 'psf_z');
psf_folder_psfmat = strcat(psf_folder,'\\','psf_zmat');
psf_folder_phase = strcat(psf_folder,'\\','phase');

if ~exist(psf_folder_H,'file')
    mkdir(psf_folder_H);
end

if ~exist(psf_folder_psf,'file')
    mkdir(psf_folder_psf);
end

if ~exist(psf_folder_psfmat,'file')
    mkdir(psf_folder_psfmat);
end

if ~exist(psf_folder_phase,'file')
    mkdir(psf_folder_phase);
end

%% 

% define object space
obj_space_x1 = 0; % point source coordination x on objective plane
obj_space_x2 = 0; % point source coordination y on objective plane
obj_space_x3 = zmin:zspacing:zmax; % point source coordination z on objective plane


obj_space = ones(length(obj_space_x1),length(obj_space_x2),length(obj_space_x3));% discrete object space

% test PSF on line on the NATIVE image plane

p3max = max(abs(obj_space_x3));
max_halfnum_ml = 28;
x1testspace = (pixelPitch/OSR)* (0:1: Nnum*OSR*max_halfnum_ml);  
[psfLine] = calcPSFFT(p3max, fobj, NA, x1testspace, lambda, M, n); %% Debays propogation
outArea = find(psfLine<Theshold);
if isempty(outArea)  
   error('Estimated PSF size exceeds the limit');   
end

half_ml_num = ceil(outArea(1)/(OSR*Nnum)); % Half of microlens is (26 + 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IMGsize = Nnum*(2*(half_ml_num+1)+1);
IMGsize_L = Nnum*(2*(half_ml_num+1)+1);



validpts = find(obj_space>eqtol); % find non-zero points
numpts = length(validpts); %
[p1indALL, p2indALL, p3indALL] = ind2sub( size(obj_space), validpts); % index to subcripts
p1ALL = obj_space_x1(p1indALL)'; % effective obj points x location
p2ALL = obj_space_x2(p2indALL)'; % effective obj points y location
p3ALL = obj_space_x3(p3indALL)'; % effective obj points z location


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% projection from points on z axis %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pixelPitch_OSR = MLPitch/OSR/Nnum; %simulated pixel size after OSR
fov=length(-(half_ml_num+1)*OSR*Nnum:(half_ml_num+1)*OSR*Nnum)*pixelPitch_OSR;   %the size of field of view for the PSF
pixelSize_OSR=length(-(half_ml_num+1)*OSR*Nnum:(half_ml_num+1)*OSR*Nnum); %the number of the pixels for each PSF

% calculate frequency domain coordinates
sinalpha_max = NA / n / M; % n sin a = NA
fx_sinalpha = 1/(2*pixelPitch_OSR);
fx_step = 1/fov ;
fx_max = fx_sinalpha ;
fx= -fx_max+fx_step/2 : fx_step : fx_max;  %% length(fx) == length(sensor)
[fx_coordination, fy_coordination] = meshgrid( fx , fx );

% calculate circle mask of optical system (Pupil range of optical system)  
fx2_coordination=fx_coordination.*fx_coordination;
fy2_coordiantion=fy_coordination.*fy_coordination;
aperture_mask=((fx2_coordination+fy2_coordiantion)<=((NA/(lambda*M)).^2));
aperture_center_line = aperture_mask((size(aperture_mask,2)+1)/2,:);
pupil_diameter = sum(aperture_center_line(:)); 

ml_space_x1 = (pixelPitch/OSR)* (-(Nnum*OSR-1)/2 : 1 : (Nnum*OSR-1)/2); % total x space per ML
ml_space_x2 = (pixelPitch/OSR)* (-(Nnum*OSR-1)/2 : 1 : (Nnum*OSR-1)/2); % total y space per ML
ml_length_x1 = length(ml_space_x1); %
ml_length_x2 = length(ml_space_x2); %

sensor_space_x1_w_osr = (pixelPitch/OSR)*(-(half_ml_num+1)*Nnum*OSR:1:(half_ml_num+1)*Nnum*OSR); % x space
sensor_space_x2_w_osr = (pixelPitch/OSR)*(-(half_ml_num+1)*Nnum*OSR:1:(half_ml_num+1)*Nnum*OSR); % y space
sensor_length_x1_osr = length(sensor_space_x1_w_osr); % x
sensor_length_x2_osr = length(sensor_space_x2_w_osr); % y

ml_array_modulation = calcML(fml, k0, ml_space_x1, ml_space_x2, sensor_space_x1_w_osr, sensor_space_x2_w_osr); % micro array phase mask modulation
ml_array_modulation = gpuArray(single(ml_array_modulation));

obj_space_x1 = (pixelPitch/M)*(-floor(Nnum/2):1:floor(Nnum/2)); % corresponding object plane coordination x1  (Nnum x Num point 
obj_space_x2 = (pixelPitch/M)*(-floor(Nnum/2):1:floor(Nnum/2)); % corresponding object plane coordination x2   within one microlens)
XREF = ceil(length(obj_space_x1)/2);
YREF = ceil(length(obj_space_x1)/2);
center_pixel_coordinate = ceil(length(sensor_space_x1_w_osr)/2);
halfWidth = (half_ml_num)*Nnum*OSR;

pixel_coordinate = ( (center_pixel_coordinate-1)/OSR+1 - halfWidth/OSR :1: (center_pixel_coordinate-1)/OSR+1 + halfWidth/OSR  );%
H_z = zeros(length(pixel_coordinate),length(pixel_coordinate), Nnum, Nnum);

psf_z = zeros(IMGsize,IMGsize,Nnum,Nnum,'single');

psf_wave_before_modulation = ones(pixelSize_OSR,pixelSize_OSR).*aperture_mask;
psf_wave_before_modulation = gpuArray(single(psf_wave_before_modulation));

%% Phase generation

if aberration_flag  == 1
    
    aberration_matrix = create_aberration(pupil_diameter,abs(ratio));
%     load(sprintf('D:\\03project\\01RUSH3D\\RUSH3Dproject\\RUSH3Dresult\\beads\\phaseiter_%d\\phase.mat',0));
%     aberration_matrix = zeros(size(phase));
%     for i = 1:iter
%         load(sprintf('D:\\03project\\01RUSH3D\\RUSH3Dproject\\RUSH3Dresult\\beads\\phaseiter_%d\\phase.mat',iter-1));
%         aberration_matrix = aberration_matrix + phase;
%     end
    save([psf_folder_phase,'\\','phase','.mat'],'aberration_matrix','-v7.3');
    
    figure(1);
    x = linspace(-1,1,pupil_diameter);
    pcolor(x,x,wrapTo2Pi(aberration_matrix)), shading interp;colorbar;title('aberration');axis equal;
    saveas(gcf,strcat(psf_folder_phase,'\\','aberration_phase','.png'));
    close(figure(1));

end


%% Calculating PSF

disp(['Start Calculating PSF...']);

% loop for different zposition
for eachpt=1:numpts
    ztime = tic;
    disp(['calcu point #', num2str(eachpt),' of ',num2str(numpts),' in total.']);
    time_s = tic;
    p1 = p1ALL(eachpt); % object point x
    p2 = p2ALL(eachpt); % object point y
    p3 = p3ALL(eachpt); % object point z
    
    % psf on center line

    wfm_modulation=k0*n*p3*realsqrt((1-(fx_coordination.*lambda./n.*M).^2-(fy_coordination.*lambda./n.*M).^2).*aperture_mask);
    wfm_modulation = gpuArray(single(wfm_modulation));
    psf_wave_original_focal_plane_fdomain = psf_wave_before_modulation.*exp(1j*wfm_modulation);
    
    
    
    if aberration_flag  == 1
        padr = round((size(psf_wave_original_focal_plane_fdomain,1) - size(aberration_matrix,1))/2);
        aberration_matrix = padarray(aberration_matrix,[padr,padr],0,'both');
        psf_wave_original_focal_plane_fdomain = psf_wave_original_focal_plane_fdomain.* exp(1j*aberration_matrix*sign(ratio));
    end
    
   
    
    psf_wave_original_focal_plane_sdomain=fftshift(ifft2(ifftshift(squeeze(psf_wave_original_focal_plane_fdomain))));
    

    timeFre = tic;
    
    % calculate for psf for each point at different (xy) position
    for b1 = 1:length(obj_space_x2)
        for a1 = 1:length(obj_space_x1)
            psf_shift_before_ml= im_shift2_GPU(psf_wave_original_focal_plane_sdomain, OSR*(a1-XREF), OSR*(b1-YREF) ); % shift positive
            psf_on_sensor_w_osr =fresnel2D_GPU(psf_shift_before_ml.*ml_array_modulation, pixelPitch/OSR, fml,lambda); % fresnel propogation for fml
            psf_on_sensor_w_osr= im_shift2_GPU(psf_on_sensor_w_osr, -OSR*(a1-XREF), -OSR*(b1-YREF) ); % shift negative (back to origin)
            
            [psf_on_sensor, x1shift, x2shift] = pixelBinning_GPU(abs(psf_on_sensor_w_osr).^2, OSR); % pixel binning
            temp_H = psf_on_sensor( pixel_coordinate - x1shift, pixel_coordinate-x2shift );
            
            H_z(:,:,a1,b1) = gather(temp_H);%
        end
    end
    
    tt = toc(timeFre);
    disp(['calcu H take ',num2str(tt),' sec....']);
    
    % remove small value
    tol_value = tol*max(H_z(:));
    H_z(H_z< (tol_value)) = 0;
    
    % normalize
    for b2 = 1:length(obj_space_x2)
        for a2 = 1:length(obj_space_x1)
            sss = H_z(:,:,a2,b2);
            H_z(:,:,a2,b2) = H_z(:,:,a2,b2)./sum(sss(:));
        end
    end
    
    %% OSR = 5 
    if OSR == 5
        temp = zeros(size(H_z,1)+1,size(H_z,2)+1,size(H_z,3),size(H_z,4));
        temp(2:end,2:end,:,:) = H_z;
        H_z = temp(1:end-1,1:end-1,:,:,:);
    end
    disp('split');
    
    border=fix(IMGsize_L/2)-fix(size(H_z,1)/2);
    blur_image=zeros(IMGsize_L,IMGsize_L,size(H_z,3),size(H_z,4));
    
    for i = 1:size(H_z,3)
        for j = 1:size(H_z,4)
            temp=zeros(IMGsize_L,IMGsize_L);
            temp(border+1:end-border,border+1:end-border)=squeeze(H_z(:,:,i,j));
            blur_image(:,:,i,j)=(im_shift3d(temp,i-((Nnum+1)/2),j-((Nnum+1)/2))); % It is used for adapt to the ReAlign function
        end
    end
    
    blur_image(isnan(blur_image)) = 0;
    maxH_z  = max(blur_image(:));
    Output=uint16(blur_image./maxH_z*65535);
    for idu=1:Nnum
        for idv = 1:Nnum
            ii = (idu-1)*Nnum+idv;
            if ii==1
                imwrite(uint16(squeeze(Output(:,:,idu,idv))),strcat(psf_folder_H, '\\', 'H_z', num2str(eachpt),'.tif'));
            else
                imwrite(uint16(squeeze(Output(:,:,idu,idv))),strcat(psf_folder_H, '\\', 'H_z',num2str(eachpt),'.tif'),'WriteMode', 'append');
            end
        end
    end
    psf_temp = strcat(psf_folder,'\\','psf_temp'); % it will be covered when next z_slice
    H_file = strcat(psf_folder_H, '\\','H_z',num2str(eachpt),'.tif');
    
    command = sprintf('ReAlign %d %s %s %d %d %d %d %d %d %s %d %d %d %d',...
        Nnum, psf_temp, H_file, 0, (IMGsize_L-1)/2, (IMGsize_L-1)/2, ((IMGsize_L/Nnum)-1)/2,...
        1, 0, '15x15.conf.sk.png',1,0, Nnum, ((IMGsize_L/Nnum)-1)/2);
    system(command);
    
    bbb = zeros(IMGsize_L,IMGsize_L,Nnum,Nnum);
    for idu = 1:Nnum
        for idv = 1:Nnum
            tmp = imread(strcat(psf_temp,'_No0.tif'),(idu-1)*Nnum+idv);
            bbb(:,:,idu,idv) = double(tmp)./65535.*maxH_z;
        end
    end
    
    tmp = bbb((IMGsize_L+1)/2-(IMGsize-1)/2:(IMGsize_L+1)/2+(IMGsize-1)/2,(IMGsize_L+1)/2-(IMGsize-1)/2:(IMGsize_L+1)/2+(IMGsize-1)/2,:,:);
    
    psf_z = single(tmp);
    temp=gather(tmp);
    Output=uint16(gather(temp)./max(bbb(:)).*65535);
    for idu=1:Nnum
        for idv = 1:Nnum
            ii = (idu-1)*Nnum+idv;
            if ii==1
                imwrite(uint16(squeeze(Output(:,:,idu,idv))),strcat(psf_folder_psf,'\\','psf_sim',num2str(eachpt),'.tif'));
            else
                imwrite(uint16(squeeze(Output(:,:,idu,idv))),strcat(psf_folder_psf,'\\','psf_sim',num2str(eachpt),'.tif'),'WriteMode', 'append');
            end
        end
    end
    
    save([psf_folder_psfmat,'\\','psf_sim',num2str(eachpt),'.mat'],'psf_z','-v7.3');
    onetime = toc(ztime);
    disp(['id_z = ',num2str(eachpt),' taketime ',num2str(onetime),' sec......']);
    
    
end
end