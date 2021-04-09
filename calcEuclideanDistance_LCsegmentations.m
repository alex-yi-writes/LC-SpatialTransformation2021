%% calculate Euclidean distance: LC segmentations
clc;clear;close all
warning off
 
% paths
path_transformed = '/path/to/the/transformed/LCsegmentations/';
IDs=[]; % subject IDs
 
% load images
% aggregated LC segmentations, transformed & binarized
LC_tf_bin = spm_read_vols(spm_vol([path_transformed 'LCmask_heatmap_binary.nii'])); LC_tf_bin(find(LC_tf_bin==0))=NaN; 
% template MNI LC mask
MNImask = spm_read_vols(spm_vol(['/path/to/the/template/LCMask/mni_icbm152/mni_icbm152_LCmetaMask_MNI05_s01f_plus50_bin.nii'])); MNImask(MNImask==0)=0;
 
% identify coordinates of each voxel
[x_MNI,y_MNI,z_MNI] = ind2sub(size(MNImask),find(MNImask~=0));
[x_tf,y_tf,z_tf] = ind2sub(size(LC_tf_bin),find(~isnan(LC_tf_bin)));
 
% coordinates in double
positions_MNI = [x_MNI,y_MNI,z_MNI];
[~,idx1] = sort(positions_MNI(:,2));
slices_MNI_mask = unique(positions_MNI(:,3)); clear idx1
 
 
%% draw and check the transformed & aggregated LC segmentations in 3D space
 
close all
 
S1 = repmat([170],numel(x_tf),1);
S2 = repmat([200],numel(x_MNI),1);
 
hFig = figure();
axh = axes('Parent', hFig);
set(gca,'FontSize',25); hold on
hold(axh, 'all');
h2 = scatter3(x_MNI,y_MNI,z_MNI,S2,'d',...
    'MarkerEdgeColor','g',...
    'MarkerFaceColor',[0.4660 0.6740 0.1880]);
h1 = scatter3(x_tf,y_tf,z_tf,S1,'o',...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor',[1 0.2 0.2],'MarkerFaceAlpha',.5);
xlabel('X','FontWeight','bold')
ylabel('Y','FontWeight','bold')
zlabel('Z','FontWeight','bold')
xlim([88 107])
ylim([91 99])
zlim([46 64])
view(axh, -33, 22);
grid(axh, 'on'); 
set(gca,'FontSize',18);
legend(axh, [h1,h2], {'Transformed Masks', 'MNI Meta LC mask'});
 
%% find centroid voxel
 
% left LC
dummybase = zeros(193,229,193); % make a blank canvas
 
tmp_MNI=[x_MNI,y_MNI,z_MNI]; % tidy up the xyz coordinates from the original mask
 
indL = tmp_MNI(:,1)<mean(x_MNI); % identify coordinates of left and right
coordsL = tmp_MNI(indL,:,:);
indR = tmp_MNI(:,1)>mean(x_MNI);
coordsR = tmp_MNI(indR,:,:);
 
leftLC=dummybase;
leftLC(sub2ind(size(leftLC),coordsL(:,1),coordsL(:,2),coordsL(:,3))) = 1;
rightLC=dummybase; 
rightLC(sub2ind(size(leftLC),coordsR(:,1),coordsR(:,2),coordsR(:,3))) = 1;
 
[x_Left,y_Left,z_Left]=ind2sub(size(leftLC),find(leftLC~=0));
[x_Right,y_Right,z_Right]=ind2sub(size(rightLC),find(rightLC~=0));
 
centroid_L_mni = [round(mean(x_Left)),round(mean(y_Left)),round(mean(z_Left))];
centroid_R_mni = [round(mean(x_Right)),round(mean(y_Right)),round(mean(z_Right))];
centroid_both_mni = [round(mean(x_MNI)),round(mean(y_MNI)),round(mean(z_MNI))];
 
%% extract centroid of single-subject masks
 
centroid_L = []; centroid_R =[];
for id = 1:length(IDs)
    
    clear tmp leftLC rightLC indL indR coordsL coordsR xi_Left yi_Left zi_Left ...
        xi_Right yi_Right zi_Right
    
    LCsegs_individual_binary{id,1} = spm_read_vols(spm_vol([path_transformed 'threshold025_' IDs{id} '_conjmask_mni.nii']));
 
    dummybase = zeros(193,229,193); % make a blank canvas
    
    [xi,yi,zi] = ind2sub(size(LCsegs_individual_binary{id,1}),find(LCsegs_individual_binary{id,1}~=0));
    tmp=[xi,yi,zi]; % tidy up the xyz coordinates from the original mask
    
    indR = tmp(:,1)>97;%mean(xi); % identify coordinates of
    coordsR = tmp(indR,:,:);
    indL = tmp(:,1)<97;%mean(xi);
    coordsL = tmp(indL,:,:);
    
    leftLC=dummybase;
    leftLC(sub2ind(size(leftLC),coordsL(:,1),coordsL(:,2),coordsL(:,3))) = 1;
    
    rightLC=dummybase;
    rightLC(sub2ind(size(leftLC),coordsR(:,1),coordsR(:,2),coordsR(:,3))) = 1;
    
    [xi_Left,yi_Left,zi_Left]=ind2sub(size(leftLC),find(leftLC~=0));
    left_indivs{id,1} = [xi_Left,yi_Left,zi_Left];
    [xi_Right,yi_Right,zi_Right]=ind2sub(size(rightLC),find(rightLC~=0));
    right_indivs{id,1} = [xi_Right,yi_Right,zi_Right];
    
    centroid_L(id,:) = [round(mean(xi_Left)),round(mean(yi_Left)),round(mean(zi_Left))];
    centroid_R(id,:) = [round(mean(xi_Right)),round(mean(yi_Right)),round(mean(zi_Right))];
    centroid_both(id,:) = [round(mean(xi)),round(mean(yi)),round(mean(zi))];
end
 
 
 
%% ED from the centroid point, slicewise analysis
 
% find centre of each slice
MNI_left = [x_Left,y_Left,z_Left];
MNI_right = [x_Right,y_Right,z_Right];
MNI_slices_L = unique(z_Left);
MNI_slices_R = unique(z_Right);
 
for slices = 1:length(MNI_slices_L)
    clear tmp zind
    % sort the left
    zind=z_Left==MNI_slices_L(slices);
    tmp = MNI_left(zind',:,:);
    if size(MNI_left(zind',:,:),1)==1
    Lcentrer_mni(slices,:) = [tmp(:,1:2) MNI_slices_L(slices)];
    else
    Lcentrer_mni(slices,:) = [mean(tmp(:,1:2)) MNI_slices_L(slices)];
    end
end
 
for slices = 1:length(MNI_slices_R)
    
    clear tmp zind
    % sort the left
    zind=z_Right==MNI_slices_R(slices);
    tmp = MNI_right(zind',:,:);
    if size(MNI_right(zind',:,:),1)==1
    Rcentrer_mni(slices,:) = [tmp(:,1:2) MNI_slices_R(slices)];
    else
    Rcentrer_mni(slices,:) = [mean(tmp(:,1:2)) MNI_slices_R(slices)];
    end
    
end
 
 
% now find centres of individual masks
Lcentre_indiv=[];Rcentre_indiv=[];
for id=1:length(IDs)
    
    clear maskL maskR slicesL slicesR sl
    
    % first, left
    maskL=right_indivs{id,1}; % there's something wrong with this left and right
    slicesL=unique(maskL(:,3));
    for sl=1:length(slicesL)
        clear tmp zind
        % sort the left
        zind=maskL(:,3)==slicesL(sl);
        tmp = maskL(zind',:,:);
        if size(tmp,1)== 1
            Lcentre_indiv{id,1}(sl,:) = [tmp(:,1:2) slicesL(sl)];
        else
            Lcentre_indiv{id,1}(sl,:) = [mean(tmp(:,1:2)) slicesL(sl)];
        end
    end
    
    
    % then, right
    maskR=left_indivs{id,1}; % there's something wrong with this left and right
    slicesR=unique(maskR(:,3));
    for sl=1:length(slicesR)
        clear tmp zind
        % sort the left
        zind=maskR(:,3)==slicesR(sl);
        tmp = maskR(zind',:,:);
        if size(tmp,1)== 1
            Rcentre_indiv{id,1}(sl,:) = [tmp(:,1:2) slicesR(sl)];
        else
            Rcentre_indiv{id,1}(sl,:) = [mean(tmp(:,1:2)) slicesR(sl)];
        end
    end
    
end
 
%% calculate Euclidian distance of centroids, slicewise
 
ED_slice_L=[]; ED_slice_R=[];
for id=1:length(IDs)
    
    % find the slices that matches
    clear Lslc_tmp A B
    Lslc_tmp = Lcentre_indiv{id,1}(:,3);
    indslc1=ismember(MNI_slices_L,Lslc_tmp,'rows');
    indslc2=ismember(Lslc_tmp,MNI_slices_L,'rows');
    
    A = Lcentrer_mni(indslc1,:);
    B = Lcentre_indiv{id,1}(indslc2,:);
    for k=1:size(B,1)
        clear v
        v    = B(k,:) - A(k,:);
        ED_slice_L{id,1}(k) = sqrt(nansum(v .^ 2));
    end
    
    % find the slices that matches
    clear Rslc_tmp A B
    Rslc_tmp = Rcentre_indiv{id,1}(:,3);
    indslc1=ismember(MNI_slices_R,Rslc_tmp,'rows');
    indslc2=ismember(Rslc_tmp, MNI_slices_R,'rows');
    
    A = Rcentrer_mni(indslc1,:);
    B = Rcentre_indiv{id,1}(indslc2,:);
    
    for k=1:size(A,1)
        clear v
        v    = B(k,:) - A(k,:);
        ED_slice_R{id,1}(k) = sqrt(nansum(v .^ 2));
    end
 
end
 
 
% additional processing
ED_slice_mean_L = cellfun(@nanmean, ED_slice_L);
ED_slice_mean_R = cellfun(@nanmean, ED_slice_R);
 
%% simple stats
 
means(1,1)=nanmean(ED_slice_mean_L);
means(1,2)=nanmean(ED_slice_mean_R);
stds(1,1)=nanstd(ED_slice_mean_L);
stds(1,2)=nanstd(ED_slice_mean_R);
medians(1,1)=nanmedian(ED_slice_mean_L);
medians(1,2)=nanmedian(ED_slice_mean_R);
 
fprintf('\n left LC mean: %1.4f, STD: %1.4f, median: %1.4f \n right LC mean: %1.4f, STD: %1.4f, median: %1.4f \n', ...
    means(1,1), stds(1,1), medians(1,1), means(1,2), stds(1,2), medians(1,2))
