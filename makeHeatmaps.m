%% Check the EPI spatial transformation with landmarks on the transformed EPI

clear;clc

% preparation

% SET PATHS
path_root = '/path/to/your/landmark/images/';
path_save = '/path/to/the/folder/where/the/heatmap/will/be/saved/';
cd(path_save)

% load IDs
load('/path/to/your/IDs.mat')

% binary or heatmap?
prompt = {'Which image will you generate?: binary or heatmap'};
dlgtitle = 'Input';
definput = {'heatmap'};
dims = [1 35];
answer=inputdlg(prompt,dlgtitle,dims,definput);
if strcmp(answer{1,1},'heatmap')
binarise=0;
elseif strcmp(answer{1,1},'binary')
binarise=1;
end

%% make

ccx = 0;ccy = 0;ccz = 0; % reset the counter
MNI = spm_read_vols(spm_vol(['/Users/alex/Dropbox/Masks/MNI/mni_icbm152_t1_tal_nlin_asym_09c.nii'])); % change here accordingly to your file name
averageplot = zeros(size(MNI)); % pre-assign the dummy variable for plotting the masks. must be the same size with the group/MNI image

for id = 1:length(ids)
    
    averplotdummy = zeros(size(MNI)); % an empty matrix that landmarks will be drawn
    
    name = ids{id};
    data = spm_read_vols(spm_vol([path_root name '/data/EPILandmark_mni.nii'])); % change here accordingly to your file name
    [x y z] = ndgrid(1:size(data,1), 1:size(data,2), 1:size(data,3));
    
    % get x/y/z level of landmarks
    ycoordz = y(find(data(:) ~= 0)); uycoordz = unique(ycoordz);
    ylevel=uycoordz;
    
    xcoordz = x(find(data(:) ~= 0)); uxcoordz = unique(xcoordz);
    xlevel=uxcoordz;
    
    zcoordz = z(find(data(:) ~= 0)); uzcoordz = unique(zcoordz);
    zlevel=uzcoordz;    
        
    ccz = ccz+1;
    for zc=1:length(zlevel)
        dumz = squeeze(data(:,:,zlevel(zc)));
        averplotdummy(:,:,zlevel(zc)) = dumz; % average points and squeeze them in
    end
    
    ccy = ccy+1;
    for yc=1:length(ylevel)
        dumy = squeeze(data(:,:,ylevel(yc)));
        averplotdummy(:,:,ylevel(yc)) = dumy; % average points and squeeze them in
    end
    
    ccx = ccx+1;
    for xc=1:length(xlevel)
        dumx = squeeze(data(:,:,xlevel(xc)));
        averplotdummy(:,:,xlevel(xc)) = dumx; % average points and squeeze them in
    end
    
    averageplot = averageplot + averplotdummy; % record the points
    clear dumx dumy dumz
    
end


%% generate the overlay image in 3D space

if binarise==0
averplot_nifti = averageplot./length(ids); 
elseif binarise==1
averplot_nifti = averageplot;
averplot_nifti(find(averplot_nifti(:)>0)) = 1;
end
hdr = spm_vol([path_root ids{1} '/data/EPILandmark_mni.nii']); % pick just any header from a file
hdr.fname = [path_save 'heatmap_landmarks.nii'];
hdr.dim = size(averplot_nifti);
hdr = rmfield(hdr,'pinfo');
hdr.nii = spm_write_vol(hdr,averplot_nifti);

%% Check the structural transformation with the transformed LC segmentations

% preparation

% SET PATHS
path_root = '/path/to/your/landmark/images/';
path_save = '/path/to/the/folder/where/the/heatmap/will/be/saved/';
cd(path_save)

% load IDs
load('/path/to/your/IDs.mat')

% binary or heatmap?
prompt = {'Which image will you generate?: binary or heatmap'};
dlgtitle = 'Input';
definput = {'heatmap'};
dims = [1 35];
answer=inputdlg(prompt,dlgtitle,dims,definput);
if strcmp(answer{1,1},'heatmap')
binarise=0;
elseif strcmp(answer{1,1},'binary')
binarise=1;
end

%% make

MNI = spm_read_vols(spm_vol(['/Users/alex/Dropbox/Masks/MNI/mni_icbm152_t1_tal_nlin_asym_09c.nii'])); % change here accordingly to your file name
averageplot = zeros(size(MNI)); % pre-assign the dummy variable for plotting the masks. must be the same size with the group/MNI image

for id = 1:length(ids)
    
    averplotdummy = zeros(size(MNI));
    name = ids{id};
    data = spm_read_vols(spm_vol([path_root name '_conjMask_NN.nii'])); % change here accordingly to your file name
    [x y z] = ndgrid(1:size(data,1), 1:size(data,2), 1:size(data,3));
    
    % get x/y/z level of landmarks
    ycoordz = y(find(data(:) ~= 0)); uycoordz = unique(ycoordz);
    ylevel=uycoordz;
    
    xcoordz = x(find(data(:) ~= 0)); uxcoordz = unique(xcoordz);
    xlevel=uxcoordz;
    
    zcoordz = z(find(data(:) ~= 0)); uzcoordz = unique(zcoordz);
    zlevel=uzcoordz;
    
    for zc=1:length(zlevel)
        dumz = squeeze(data(:,:,zlevel(zc)));
        averplotdummy(:,:,zlevel(zc)) = dumz; % average points and squeeze them in
    end
    
    for yc=1:length(ylevel)
        dumy = squeeze(data(:,:,ylevel(yc)));
        averplotdummy(:,:,ylevel(yc)) = dumy; % average points and squeeze them in
    end
    
    for xc=1:length(xlevel)
        dumx = squeeze(data(:,:,xlevel(xc)));
        averplotdummy(:,:,xlevel(xc)) = dumx; % average points and squeeze them in
    end
    
    averplotdummy(find(averplotdummy(:)>0)) = 1;
    
    averageplot = averageplot + averplotdummy; % record the points

    clear dumx dumy dumz
    
end

%% generate the overlay image in 3D space

if binarise==0
averplot_nifti = averageplot./length(ids); 
elseif binarise==1
averplot_nifti = averageplot;
averplot_nifti(find(averplot_nifti(:)>0)) = 1;
end
hdr = spm_vol([path_root ids{1} '_conjMask_NN.nii']); % pick just any header from a file
hdr.fname = [path_save 'heatmap_segmentations.nii'];
hdr.dim = size(averplot_nifti);
hdr = rmfield(hdr,'pinfo');
hdr.nii = spm_write_vol(hdr,averplot_nifti);
