function [pupilPlaneMatrix_overAll,maxIncreament] = phaseReconstructionAccordingToGradients(PSFParameters, gradients_xDirection, gradients_yDirection, folderName_phaseReconstruction)
addpath(genpath('utils'));
xysize = PSFParameters.xysize;
centrePixelIndex = PSFParameters.centrePixelIndex;
kMaxPixels = PSFParameters.kMaxPixels;
shift_x_normalized = PSFParameters.shift_x_normalized;
shift_y_normalized = PSFParameters.shift_y_normalized;
previousPupilPlaneMatrix = PSFParameters.pupilPlaneMatrix;


%%
xn = linspace(centrePixelIndex-xysize,xysize-centrePixelIndex,xysize)/kMaxPixels;
yn = linspace(centrePixelIndex-xysize,xysize-centrePixelIndex,xysize)/kMaxPixels;
[xcoordinates,ycoordinates] = meshgrid(xn,yn);
pupil=Circle([xysize,xysize],[centrePixelIndex,centrePixelIndex],[0,kMaxPixels]);

Sx = griddata(shift_x_normalized, shift_y_normalized, gradients_xDirection,xcoordinates,ycoordinates,'v4');
Sy = griddata(shift_x_normalized, shift_y_normalized, gradients_yDirection,xcoordinates,ycoordinates,'v4');
Sx(isnan(Sx)) = 0;
Sy(isnan(Sy)) = 0;
Sx = Sx.*pupil;
Sy = Sy.*pupil;

pupilPlaneMatrix_supplement = hfli2(Sx, Sy, xcoordinates, ycoordinates);
pupilPlaneMatrix_supplement(isnan(pupilPlaneMatrix_supplement)) = 0;



meanPhase=mean(pupilPlaneMatrix_supplement(:));
pupilPlaneMatrix_supplement(pupil==0) = meanPhase;
pupilPlaneMatrix_supplement=pupilPlaneMatrix_supplement - meanPhase;



gcf = figure('color',[1 1 1]);
imagesc(xn,yn,pupilPlaneMatrix_supplement);
colormap gray
axis xy
axis square
%set(gca, 'TickDir', 'out')
title(['Reconstructed phase -- supplement'],'FontSize', 10);
xlabel('Normalized x pupil coordinate');
ylabel('Normalized y pupil coordinate');
%colormap gray
colormap jet
colorbar
saveas(gcf,[folderName_phaseReconstruction,'/pupilPlaneMatrix_supplement.jpg']);
save([folderName_phaseReconstruction,'/pupilPlaneMatrix_supplement.mat'],'pupilPlaneMatrix_supplement','-v7.3');


maxIncreament = max(abs(pupilPlaneMatrix_supplement(:)));
pupilPlaneMatrix_overAll = pupilPlaneMatrix_supplement*3+previousPupilPlaneMatrix;
gcf = figure('color',[1 1 1]);
imagesc(xn,yn,pupilPlaneMatrix_overAll);
colormap gray
axis xy
axis square
%set(gca, 'TickDir', 'out')
title(['Reconstructed phase -- overall'],'FontSize', 10);
xlabel('Normalized x pupil coordinate');
ylabel('Normalized y pupil coordinate');
%colormap gray
colormap jet
colorbar
saveas(gcf,[folderName_phaseReconstruction,'/pupilPlaneMatrix_overAll.jpg']);
save([folderName_phaseReconstruction,'/pupilPlaneMatrix_overAll.mat'],'pupilPlaneMatrix_overAll','-v7.3');


close all;
end