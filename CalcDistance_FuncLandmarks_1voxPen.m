%% Calculate the in-plane distances


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% use this script when using 1-voxel seg pen instead of 3-voxel spherical pen 
%                            ^^^^^^^
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;clear;close all

% who rated the images?
raterLabel=inputdlg('Who evaluated the landmarks? Type initials without spaces (e.g. AY)');

% set env
path_transformed = '/path/to/your/transformed/mean/functional/MRI/landmarks/';
load('/load/list/of/IDs/in/cell/structure/ID.mat')

% load individual landmark images
transformed_landmarks_single=[];
for subj=1:length(IDs)
    transformed_landmarks_single{subj,1}=spm_read_vols(spm_vol([ path_transformed IDs{subj} '_landmark.nii']));
end; clear subj
 
% load the landmark image drawn on MNI
MNIlandmark = spm_read_vols(spm_vol(['/path/to/your/MNI-landmarks/MNI_landmarks_v7_2Labels.nii'])); % change here accordingly to your file name
MNIlandmark(find(MNIlandmark<2))=NaN;
 
% write down the coordinates of MNI landmarks
[x_mni,y_mni,z_mni] = ind2sub(size(MNIlandmark),find(~isnan(MNIlandmark)));
tmp = [x_mni,y_mni,z_mni];
[~,idx1] = sort(tmp(:,3)); tmp_sorted = tmp(idx1,:);

% sort the coordinates by slices
slices = unique(tmp(:,3));
slices = slices(end:-1:1);

% the first slice: periaqueductal grey
clear idx1 positionSlice positionLR_sorted cluster1 cluster2 cluster3 TopSliceMarker ordersLR

positionSlice = tmp_sorted(find(tmp_sorted(:,3)==slices(1)),:); % find the points in the first slice
[~,idx1] = sort(positionSlice(:,1)); positionLR_sorted = positionSlice(idx1,:);
TopSliceMarker = positionLR_sorted;

% periaqueductal grey
MNI.PeriaqueductalGrey = TopSliceMarker; %ceil( median(cluster3) );

% the second slice, with four landmarks: brainstem outline and 4th
% ventricle borders

clear idx1 positionSlice positionLR_sorted cluster1 cluster2 cluster3 MiddleSliceClusters ordersLR

positionSlice = tmp_sorted(find(tmp_sorted(:,3)==slices(2)),:);
[~,idx1] = sort(positionSlice(:,1)); positionLR_sorted = positionSlice(idx1,:);
MiddleSliceClusters = kmeans(positionLR_sorted(:,1),4);
ordersLR=unique(MiddleSliceClusters,'stable'); % identify cluster indice

% brainstem outline right
cluster1 = positionLR_sorted(MiddleSliceClusters==ordersLR(4),:);
MNI.OutlineBrainstem_right = [median(cluster1(:,1)) median(cluster1(:,2)) slices(2)]; %ceil( median(cluster1) );

% 4th ventricle border right (potential LC location)
cluster2 = positionLR_sorted(MiddleSliceClusters==ordersLR(3),:);
MNI.LC_right = [median(cluster2(:,1)) median(cluster2(:,2)) slices(2)];

% 4th ventricle border left (potential LC location)
cluster3 = positionLR_sorted(MiddleSliceClusters==ordersLR(2),:);
MNI.LC_left = [median(cluster3(:,1)) median(cluster3(:,2)) slices(2)];

% brainstem outline left
cluster4 = positionLR_sorted(MiddleSliceClusters==ordersLR(1),:);
MNI.OutlineBrainstem_left = [median(cluster4(:,1)) median(cluster4(:,2)) slices(2)];


% the third slice, with one landmark: bottom window brainstem (perifastigial sulcus)

clear idx1 positionSlice positionLR_sorted cluster1 cluster2 cluster3 BottomSliceMarker

positionSlice = tmp_sorted(find(tmp_sorted(:,3)==slices(3)),:);
[~,idx1] = sort(positionSlice(:,1)); positionLR_sorted = positionSlice(idx1,:);
BottomSliceMarker = positionLR_sorted;

% perifastigial sulcus
MNI.PerifastigialSulcus = BottomSliceMarker;

disp('prep done')
 
%% identify median coords in the single landmarks
 
transformed_landmark_coords=[];
 
for subj=1:length(IDs)
    
    clear xs ys zs tmp idx2
    
    [xs,ys,zs] = ind2sub(size(transformed_landmarks_single{subj}),find(transformed_landmarks_single{subj}~=0)); % index the coordinates
    tmp = [xs ys zs]; % 2-dim matrix
    [~,idx2] = sort(tmp(:,3)); tmp_sorted = tmp(idx2,:);
    
    slices = unique(tmp(:,3)); % sort with slice numbers
    slices = slices(end:-1:1); % order into S -> I direction
    
    %% the first slice: periaqueductal grey
    
    clear idx3 positionSlice positionX_sorted
    
    transformed_landmark_coords{subj,1}.TopSlice = [];
    
    positionSlice = tmp_sorted(find(tmp_sorted(:,3)==slices(1)),:); % find the points in the first slice
    [~,idx3] = sort(positionSlice(:,1)); positionLR_sorted = positionSlice(idx3,:); % sort coordinates in L -> R direction
        
    % periaqueductal grey
    clear r c
    [r,c]=find(positionLR_sorted(:,1)==median(positionLR_sorted(:,1),'all')); % find the voxel that are in the middle
    transformed_landmark_coords{subj,1}.TopSlice.PeriaqueductalGrey = positionLR_sorted(r,:);
    
    %% the second slice, with four landmarks: brainstem outline and 4th ventricle borders
    
    clear idx3 idx4 positionY positionX_sorted
    
    transformed_landmark_coords{subj,1}.MidSlice = [];
    
    positionSlice = tmp_sorted(find(tmp_sorted(:,3)==slices(2)),:); % find the points in the second slice
    [~,idx3] = sort(positionSlice(:,1)); positionLR_sorted = positionSlice(idx3,:); % sort coordinates in L -> R direction
    
    % brainstem outline left
    clear r c
    [r,c]=find(positionLR_sorted(:,1)==positionLR_sorted(1,1));
    transformed_landmark_coords{subj,1}.MidSlice.OutlineBrainstem_left = positionLR_sorted(r,:);

    % 4th ventricle border right
    clear r c
    [r,c]=find(positionLR_sorted(:,1)==positionLR_sorted(2,1));
    transformed_landmark_coords{subj,1}.MidSlice.LC_left = positionLR_sorted(r,:);
    
    % 4th ventricle border left
    clear r c
    [r,c]=find(positionLR_sorted(:,1)==positionLR_sorted(3,1));
    transformed_landmark_coords{subj,1}.MidSlice.LC_right = positionLR_sorted(r,:);
    
    % brainstem outline right
    clear r c
    [r,c]=find(positionLR_sorted(:,1)==positionLR_sorted(4,1));
    transformed_landmark_coords{subj,1}.MidSlice.OutlineBrainstem_right = positionLR_sorted(r,:);
    
    
    %% the third slice, with one landmark: bottom window brainstem
    
    clear idx3 positionY positionX_sorted
    
    transformed_landmark_coords{subj,1}.BottomSlice = [];
    
    positionSlice = tmp_sorted(find(tmp_sorted(:,3)==slices(3)),:); % find the points in the third slice
    [~,idx3] = sort(positionSlice(:,1)); positionLR_sorted = positionSlice(idx3,:);
    
    % perifastigial sulcus
    transformed_landmark_coords{subj,1}.BottomSlice.PerifastigialSulcus = positionLR_sorted;
    
    fprintf('\n subject %s done\n',IDs{subj})
    
end;clear subj
 
disp('coordinates collected')
 
%% calculate distances in single subject level

% flexible listing
varnames = {'PeriaqueductalGrey';'OutlineBrainstem_left';'OutlineBrainstem_right';...
    'LC_left';'LC_right';'PerifastigialSulcus'};
for v1=1:length(varnames)
   eval([varnames{v1} '=[];']) 
end
 
Distances_indv=[]; 
for subj=1:length(IDs)
    
    Distances_indv{subj,1}.PeriaqueductalGrey = sum((transformed_landmark_coords{subj,1}.TopSlice.PeriaqueductalGrey-MNI.PeriaqueductalGrey).^2).^0.5;
    
    Distances_indv{subj,1}.OutlineBrainstem=[];
    Distances_indv{subj,1}.OutlineBrainstem.left=sum((transformed_landmark_coords{subj,1}.MidSlice.OutlineBrainstem_left-MNI.OutlineBrainstem_left).^2).^0.5;
    Distances_indv{subj,1}.OutlineBrainstem.right=sum((transformed_landmark_coords{subj,1}.MidSlice.OutlineBrainstem_right-MNI.OutlineBrainstem_right).^2).^0.5;
    
    Distances_indv{subj,1}.LC=[];
    Distances_indv{subj,1}.LC.left=sum((transformed_landmark_coords{subj,1}.MidSlice.LC_left-MNI.LC_left).^2).^0.5;
    Distances_indv{subj,1}.LC.right=sum((transformed_landmark_coords{subj,1}.MidSlice.LC_right-MNI.LC_right).^2).^0.5;
    
    Distances_indv{subj,1}.PerifastigialSulcus = sum((transformed_landmark_coords{subj,1}.BottomSlice.PerifastigialSulcus-MNI.PerifastigialSulcus).^2).^0.5;
    
    % write as a table
    PeriaqueductalGrey(subj,1)=Distances_indv{subj,1}.PeriaqueductalGrey;
    OutlineBrainstem_left(subj,1)=Distances_indv{subj,1}.OutlineBrainstem.left;
    OutlineBrainstem_right(subj,1)=Distances_indv{subj,1}.OutlineBrainstem.right;
    LC_left(subj,1)=Distances_indv{subj,1}.LC.left;
    LC_right(subj,1)=Distances_indv{subj,1}.LC.right;
    PerifastigialSulcus(subj,1)=Distances_indv{subj,1}.PerifastigialSulcus;
    
end
 
IDcolumn = cell2mat(cellfun(@str2num, IDs,'UniformOutput',false));
 
T=table(PeriaqueductalGrey,OutlineBrainstem_left,OutlineBrainstem_right,...
    LC_left,LC_right,PerifastigialSulcus);
Distance_export =[IDcolumn',PeriaqueductalGrey,OutlineBrainstem_left,OutlineBrainstem_right,...
    LC_left,LC_right,PerifastigialSulcus];
 
save([path_transformed 'Distance_EPILandmarks_' raterLabel '.mat'],'Distance_export')
save([path_transformed 'LandmarkCoordinates_' raterLabel '.mat'],'transformed_landmark_coords')
 
disp('distance calc done')

%% calculate interrator agreement

load('/load/list/of/IDs/in/cell/structure/ID.mat') % load IDs of rated images

clear transformed_landmark_coords
rater1=load('/load/the/coordinates/of/landmarks/LandmarkCoordinates_rater1.mat');
rater2=load('/load/the/coordinates/of/landmarks/LandmarkCoordinates_rater2.mat');

% calculate agreements per subject per landmark
IR=[];
for id=1:length(IDs)
    
        
    % top slice
    clear rater1dat rater2dat
    rater1dat=rater1.transformed_landmark_coords{id,1}.TopSlice;
    rater2dat=rater2.transformed_landmark_coords{id,1}.TopSlice;
    
    IR{id,1}.PeriaqueductalGrey       =sum(rater1dat.PeriaqueductalGrey(1:2)==rater2dat.PeriaqueductalGrey(1:2));
    PeriaqueductalGrey(id,1)          =sum(rater1dat.PeriaqueductalGrey(1:2)==rater2dat.PeriaqueductalGrey(1:2));
    
    cRater1.PeriaqueductalGrey(id,:)          = rater1dat.PeriaqueductalGrey;
    cRater2.PeriaqueductalGrey(id,:)          = rater2dat.PeriaqueductalGrey;
    
    
    % middle slice
    clear ay ml
    rater1dat=rater1.transformed_landmark_coords{id,1}.MidSlice;
    rater2dat=rater2.transformed_landmark_coords{id,1}.MidSlice;
    IR{id,1}.OutlineBrainstem_left      =sum(rater1dat.OutlineBrainstem_left(1:2)==rater2dat.OutlineBrainstem_left(1:2));
    IR{id,1}.OutlineBrainstem_right     =sum(rater1dat.OutlineBrainstem_right(1:2)==rater2dat.OutlineBrainstem_right(1:2));
    IR{id,1}.LC_left                    =sum(rater1dat.LC_left(1:2)==rater2dat.LC_left(1:2));
    IR{id,1}.LC_right                   =sum(rater1dat.LC_right(1:2)==rater2dat.LC_right(1:2));
    
    OutlineBrainstem_left(id,1)       =sum(rater1dat.OutlineBrainstem_left(1:2)==rater2dat.OutlineBrainstem_left(1:2));
    OutlineBrainstem_right(id,1)      =sum(rater1dat.OutlineBrainstem_right(1:2)==rater2dat.OutlineBrainstem_right(1:2));
    LC_left(id,1)                     =sum(rater1dat.LC_left(1:2)==rater2dat.LC_left(1:2));
    LC_right(id,1)                    =sum(rater1dat.LC_right(1:2)==rater2dat.LC_right(1:2));
    
    
    cRater1.OutlineBrainstem_left(id,:) = rater1dat.OutlineBrainstem_left;
    cRater1.OutlineBrainstem_right(id,:)= rater1dat.OutlineBrainstem_right;
    cRater1.LC_left(id,:)               = rater1dat.LC_left;
    cRater1.LC_right(id,:)              = rater1dat.LC_right;
        
    cRater2.OutlineBrainstem_left(id,:) = rater2dat.OutlineBrainstem_left;
    cRater2.OutlineBrainstem_right(id,:)= rater2dat.OutlineBrainstem_right;
    cRater2.LC_left(id,:)               = rater2dat.LC_left;
    cRater2.LC_right(id,:)              = rater2dat.LC_right;
    
    
    % bottom slice
    clear ay ml
    rater1dat=rater1.transformed_landmark_coords{id,1}.BottomSlice;
    rater2dat=rater2.transformed_landmark_coords{id,1}.BottomSlice;
    IR{id,1}.PerifastigialSulcus        =sum(rater1dat.PerifastigialSulcus(1:2)==rater2dat.PerifastigialSulcus(1:2));
    PerifastigialSulcus(id,1)           =sum(rater1dat.PerifastigialSulcus(1:2)==rater2dat.PerifastigialSulcus(1:2));

    
    cRater1.PerifastigialSulcus(id,:)       = rater1dat.PerifastigialSulcus;
    cRater2.PerifastigialSulcus(id,:)       = rater2dat.PerifastigialSulcus;
    
end

% calculate DICE score

DICE=[];

DICE.PeriaqueductalGrey=(sum(PeriaqueductalGrey)*2)/(2*2*numel(IDs));

DICE.OutlineBrainstem_left=(sum(OutlineBrainstem_left)*2)/(2*2*numel(IDs));
DICE.OutlineBrainstem_right=(sum(OutlineBrainstem_right)*2)/(2*2*numel(IDs));
DICE.LC_left=(sum(LC_left)*2)/(2*2*numel(IDs));
DICE.LC_right=(sum(LC_right)*2)/(2*2*numel(IDs));

DICE.PerifastigialSulcus=(sum(PerifastigialSulcus)*2)/(2*2*numel(IDs));
