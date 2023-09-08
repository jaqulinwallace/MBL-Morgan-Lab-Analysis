% -------------------------------------------------------------------------
% Compare intensities from different ROIs
% Morgan lab [JW 2023]
% -------------------------------------------------------------------------
clear all
clc
% -------------------------------------------------------------------------
%  User Variables:
Ir = imread(''); % Red channel (594)
Ig = imread(''); % Green channel (488)
rad = 4; % radius for disk erosion of synapse edges (pixels)
minarea = 150; % minimum area required to be considered an object of interest (pixels)
include_background = 0; % Include a background subtraction 1=yes 0=no
q = 23.5211; % pixel per micron ratio
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Code starts here:
% Ir = double(Ir);
Ig = double(Ig);

if include_background == 1
    % background mask:
    imshow(Ig)
    title('Draw Background Mask - Rectangle')
    imcontrast
    mask_rect = drawrectangle;
    Ig_bckgrd_mask = createMask(mask_rect);
    close
    Ig_bckgrd_mask = double(Ig_bckgrd_mask);
    Ig_bckgrd_mask(Ig_bckgrd_mask==0) = -1;
    Ig_bckgrd = Ig.*Ig_bckgrd_mask;
    Ig_bckgrd(Ig_bckgrd<0) = [];
    bckgrd = ceil(mean(Ig_bckgrd(Ig_bckgrd>0)));
    
    % Background subtraction:
    Ig(Ig<=bckgrd) = 0;
end

% Draw axonal mask (from green channel):
imshow(Ig)
title('Draw Axonal Mask - Freehand')
imcontrast
mask_axon = drawfreehand;
Ig_axon_mask = createMask(mask_axon); % Axonal mask
close

% Threshold the red channel to isolate synapses:
kmask = imsegkmeans(Ir,2);
% Order by clusters by signal strength
if mean(Ir(kmask==1)) > mean(Ir(kmask==2))
    kmask = kmask==1;
else
    kmask = kmask==2;
end
kmask = bwareaopen(kmask,minarea);
kmask = imfill(kmask,'holes');
Ir_syn_mask = imopen(kmask,strel('disk',rad)); % Synapses mask

% Isolate synapses inside Axon and make mask:
Ir_axon = Ir_syn_mask.*Ig_axon_mask;
Ir_syn_axon_mask = Ir_axon;
Ir_syn_axon_mask(Ir_syn_axon_mask>0) = 1;  % Synapses inside Axonal Mask 

% Isolate axoplasm inside axon and make mask:
Ir_axoplasm_mask = ~Ir_syn_axon_mask.*Ig_axon_mask;

% Intensity inside axon (intensity counts):
Ig_axon = Ig.*Ig_axon_mask;
Axon_intensity = sum(Ig_axon(:));
% Intensity inside synapses (intensity counts):
Ig_axon_syn = Ig.*Ir_syn_axon_mask;
Synapse_intensity = sum(Ig_axon_syn(:))
% Intensity inside axoplasm (intensity counts):
Ig_axoplasm = Ig.*Ir_axoplasm_mask;
Axoplasm_intensity = sum(Ig_axoplasm(:))

q2 = q*q; % area in micron squared of each pixel
% Area inside axon (microns):
Axon_area = sum(Ig_axon_mask(:))./q2;
% Area inside synapses (microns):
Synapse_area = sum(Ir_syn_axon_mask(:))./q2;
% Area inside axoplasm (microns):
Axoplasm_area = sum(Ir_axoplasm_mask(:))./q2;

% Intensity density inside synapses:
Synapse_density = Synapse_intensity/Synapse_area;
% Intensity density inside axoplasm:
Axoplasm_density = Axoplasm_intensity/Axoplasm_area;

% Save workspace:
timestamp=datestr(now,'mm-dd-yy+HH-MM-SS');
save(['',timestamp,'.mat'])
% -------------------------------------------------------------------------