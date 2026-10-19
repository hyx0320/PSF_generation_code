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

M =         5.761; % M = NA*2*fml/MLPitch;% magnification || formula(15) ?5.6509
Nnum =      15; % 15 virtual pixels per microlens
n =         1; % refractive index
NA =        0.5; % NA of objective lens is 0.5.
MLPitch =   69*1e-6; % The pitch size of microlens is 69 um
fml =       393.3*1e-6; % microlens focal length 193 um
ftl =       265e-3; % Focal length of tube lens is 165 mm
OSR =       5; % calculate 5 times per pixel
lambda =    525*1e-9; %  light wave length is 525 nm

% -------------- Axial range and sampling rate ---------------- %

zmax_u = 200;
zmin_u = -200;
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

%% Save PSF path

psf_folder1 = 'D:\RUSH3Dproject\RUSH3Dpsf'; 
psf_folder_parameter = strcat('WF',num2str(M),'x_lambda',num2str(lambda*1e9),'_axial_',num2str(zmin_u),'um_to_',num2str(zmax_u), ...
                    'um_zstep_',num2str(zspacing*1e6),'_NA_',num2str(NA),'_fml_', ...
                    num2str(fml*1e6),'um_ftl_',num2str(ftl*1e3),'mm_OSR_',num2str(OSR));

psf_folder = strcat(psf_folder1, '\\', psf_folder_parameter);
psf_folder_H = strcat(psf_folder,'\\', 'H_z');
psf_folder_psf = strcat(psf_folder,'\\', 'psf_z');
psf_folder_psfmat = strcat(psf_folder,'\\','psf_zmat');


if ~exist(psf_folder_H,'file')
    mkdir(psf_folder_H);
end

if ~exist(psf_folder_psf,'file')
    mkdir(psf_folder_psf);
end

if ~exist(psf_folder_psfmat,'file')
    mkdir(psf_folder_psfmat);
end

%% 

% define object space
obj_space_x1 = 0; % point source coordination x on objective plane
obj_space_x2 = 0; % point source coordination y on objective plane
obj_space_x3 = zmin:zspacing:zmax; % point source coordination z on objective plane


obj_space = ones(length(obj_space_x1),length(obj_space_x2),length(obj_space_x3));% discrete object space

% test PSF on line on the NATIVE image plane

p3max = max(abs(obj_space_x3));
max_halfnum_ml = 85;
x1testspace = (pixelPitch/OSR)* (0:1: Nnum*OSR*max_halfnum_ml);  
[psfLine] = calcPSFFT(p3max, fobj, NA, x1testspace, lambda, M, n); %% Debays propogation
outArea = find(psfLine<Theshold);
if isempty(outArea)   %%¹ÀËã´ó¸ÅÑÓÉêµ½¶àÔ¶
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
H_z = zeros(length(pixel_coordinate),length(pixel_coordinate),Nnum, Nnum);
WF_z = zeros(length(pixel_coordinate),length(pixel_coordinate));

psf_z = zeros(IMGsize,IMGsize,Nnum,Nnum,'single');

psf_wave_before_modulation = ones(pixelSize_OSR,pixelSize_OSR).*aperture_mask;
psf_wave_before_modulation = gpuArray(single(psf_wave_before_modulation));




%% Calculating PSF

disp(['Start Calculating PSF...']);

% loop for different z-position
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
    psf_wave_original_focal_plane_sdomain=fftshift(ifft2(ifftshift(squeeze(psf_wave_original_focal_plane_fdomain))));
    

    timeFre = tic;
    
    % calculate for psf for each point at different (xy) position

    [psf_on_sensor, x1shift, x2shift] = pixelBinning_GPU(abs(psf_wave_original_focal_plane_sdomain).^2, OSR); % pixel binning
    temp_H = psf_on_sensor( pixel_coordinate - x1shift, pixel_coordinate-x2shift );
    
    WF_z = gather(temp_H);%

    tt = toc(timeFre);
    disp(['calcu H take ',num2str(tt),' sec....']);
    
    % remove small value
    tol_value = tol*max(WF_z(:));
    WF_z(WF_z< (tol_value)) = 0;
    
    % normalize
    
    WF_z = WF_z./sum(WF_z(:));

    %% if OSR == 5 
    if OSR == 5
        temp = zeros(size(WF_z,1)+1,size(WF_z,2)+1);
        temp(2:end,2:end) = WF_z;
        WF_z = temp(1:end-1,1:end-1);
    end
    disp('split');
    
    border=fix(IMGsize_L/2)-fix(size(WF_z,1)/2);
%     blur_image=zeros(IMGsize_L,IMGsize_L);
    
    
    temp=zeros(IMGsize_L,IMGsize_L);
    temp(border+1:end-border,border+1:end-border)=squeeze(WF_z(:,:));
    blur_image = temp; % It is used for adapt to the ReAlign function

    
    blur_image(isnan(blur_image)) = 0;
    Output=im2uint16(blur_image.*10);

    imwrite(uint16(squeeze(Output)),strcat(psf_folder_H, '\\', 'WF_z', num2str(eachpt),'.tif'));
    
    save([psf_folder_psfmat,'\\','psf_sim',num2str(eachpt),'.mat'],'blur_image','-v7.3');
    onetime = toc(ztime);
    disp(['id_z = ',num2str(eachpt),' taketime ',num2str(onetime),' sec......']);
    
    
end