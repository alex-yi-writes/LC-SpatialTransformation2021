%% Calculate the in-plane distances for transformed mean functional image landmarks


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% use this script when using size-3-spherical seg pen instead of 1-voxel pen 
%                            ^^^^^^^^^^^^^^^^
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;clear;close all
 
% load IDs and set a path to work on your data
path_transformed = '/path/to/your/data/';
load('/path/to/your/subject/ID/list/ID.mat')

% load individual landmark images
transformed_landmarks_single=[];
for subj=1:length(IDs)
    transformed_landmarks_single{subj,1}=spm_read_vols(spm_vol([ path_transformed IDs{subj} '_landmark.nii']));
end; clear subj
 

%% Define MNI-based landmarks

MNI=[];

% load the landmark image drawn on MNI
MNIlandmark = spm_read_vols(spm_vol(['/path/to/your/MNI-landmarks/MNI_landmarks_for_EPI.nii'])); % change here accordingly to your file name
MNIlandmark(find(MNIlandmark<1))=NaN;
 
% write down the coordinates of MNI landmarks
[x_mni,y_mni,z_mni] = ind2sub(size(MNIlandmark),find(~isnan(MNIlandmark)));
tmp = [x_mni,y_mni,z_mni];
[~,idx1] = sort(tmp(:,3)); tmp_sorted = tmp(idx1,:);

% sort the coordinates by slices
slices = unique(tmp(:,3));
slices = slices(end:-1:1);

% the first slice, with three landmarks: top, NucRubs
clear idx1 positionSlice positionLR_sorted cluster1 cluster2 cluster3 TopSliceClusters ordersLR

positionSlice = tmp_sorted(find(tmp_sorted(:,3)==slices(1)),:); % find the points in the first slice
[~,idx1] = sort(positionSlice(:,1)); positionLR_sorted = positionSlice(idx1,:);
TopSliceClusters = kmeans(positionLR_sorted(:,1),3);
ordersLR=unique(TopSliceClusters,'stable');

% ruber right
cluster1 = positionLR_sorted(TopSliceClusters==ordersLR(3),:);
MNI.NucRuber_right = [median(cluster1(:,1)) median(cluster1(:,2)) slices(1)]; %ceil( median(cluster3) );

% top curve
cluster2 = positionLR_sorted(TopSliceClusters==ordersLR(2),:);
MNI.TopBrainstem = [median(cluster2(:,1)) median(cluster2(:,2)) slices(1)]; %ceil( median(cluster3) );

% ruber left
cluster3 = positionLR_sorted(TopSliceClusters==ordersLR(1),:);
MNI.NucRuber_left = [median(cluster3(:,1)) median(cluster3(:,2)) slices(1)]; %ceil( median(cluster3) );

% the second slice, with four landmarks: brainstem outline and LC

clear idx1 positionSlice positionLR_sorted cluster1 cluster2 cluster3 MiddleSliceClusters ordersLR

positionSlice = tmp_sorted(find(tmp_sorted(:,3)==slices(2)),:);
[~,idx1] = sort(positionSlice(:,1)); positionLR_sorted = positionSlice(idx1,:);
MiddleSliceClusters = kmeans(positionLR_sorted(:,1),4);
ordersLR=unique(MiddleSliceClusters,'stable'); % identify cluster indice

% brainstem outline right
cluster1 = positionLR_sorted(MiddleSliceClusters==ordersLR(4),:);
MNI.OutlineBrainstem_right = [median(cluster1(:,1)) median(cluster1(:,2)) slices(2)]; %ceil( median(cluster1) );

% LC right
cluster2 = positionLR_sorted(MiddleSliceClusters==ordersLR(3),:);
MNI.LC_right = [median(cluster2(:,1)) median(cluster2(:,2)) slices(2)];

% LC left
cluster3 = positionLR_sorted(MiddleSliceClusters==ordersLR(2),:);
MNI.LC_left = [median(cluster3(:,1)) median(cluster3(:,2)) slices(2)];

% brainstem outline left
cluster4 = positionLR_sorted(MiddleSliceClusters==ordersLR(1),:);
MNI.OutlineBrainstem_left = [median(cluster4(:,1)) median(cluster4(:,2)) slices(2)];


% the third slice, with one landmark: bottom window brainstem

clear idx1 positionSlice positionLR_sorted cluster1 cluster2 cluster3 BottomSliceClusters

positionSlice = tmp_sorted(find(tmp_sorted(:,3)==slices(3)),:);
[~,idx1] = sort(positionSlice(:,1)); positionLR_sorted = positionSlice(idx1,:);

% brainstem bottom
MNI.BottomBrainstem = [median(positionLR_sorted(:,1)) median(positionLR_sorted(:,2)) slices(3)]; %ceil( median(positionLR_sorted) );
 
disp('prep done')
 
%% identify median coords in the single landmarks
 
transformed_landmark_coords=[];
 
for subj=1:length(IDs)
    
    clear xs ys zs tmp idx2
    
    [xs,ys,zs] = ind2sub(size(transformed_landmarks_single{subj}),find(transformed_landmarks_single{subj}~=0));
    tmp = [xs ys zs];
    [~,idx2] = sort(tmp(:,3)); tmp_sorted = tmp(idx2,:);
    
    slices = unique(tmp(:,3)); 
    slices = slices(end:-1:1);
    
    %% the first slice, with three landmarks: top, NucRubs
    
    clear idx3 positionSlice positionLR_sorted cluster1 cluster2 cluster3 TopSliceClusters ordersLR
    
    transformed_landmark_coords{subj,1}.TopSlice = [];
    
    positionSlice = tmp_sorted(find(tmp_sorted(:,3)==slices(1)),:); % find the points in the first slice
    [~,idx3] = sort(positionSlice(:,1)); positionLR_sorted = positionSlice(idx3,:);
    TopSliceClusters = kmeans(positionLR_sorted(:,1),3);
    ordersLR=unique(TopSliceClusters,'stable');
    
    % ruber right
    cluster1 = positionLR_sorted(TopSliceClusters==ordersLR(3),:);
    transformed_landmark_coords{subj,1}.TopSlice.NucRuber_right = [median(cluster1(:,1)) median(cluster1(:,2)) slices(1)]; %ceil( median(cluster3) );
    
    % top curve
    cluster2 = positionLR_sorted(TopSliceClusters==ordersLR(2),:);
    transformed_landmark_coords{subj,1}.TopSlice.TopBrainstem = [median(cluster2(:,1)) median(cluster2(:,2)) slices(1)]; %ceil( median(cluster3) );
    
    % ruber left
    cluster3 = positionLR_sorted(TopSliceClusters==ordersLR(1),:);
    transformed_landmark_coords{subj,1}.TopSlice.NucRuber_left = [median(cluster3(:,1)) median(cluster3(:,2)) slices(1)]; %ceil( median(cluster3) );
    
    %% the second slice, with four landmarks: brainstem outline and LC
    
    clear idx3 positionSlice positionLR_sorted cluster1 cluster2 cluster3 MiddleSliceClusters ordersLR
    
    transformed_landmark_coords{subj,1}.MidSlice = [];
    
    positionSlice = tmp_sorted(find(tmp_sorted(:,3)==slices(2)),:);
    [~,idx3] = sort(positionSlice(:,1)); positionLR_sorted = positionSlice(idx3,:);
    MiddleSliceClusters = kmeans(positionLR_sorted(:,1),4);
    ordersLR=unique(MiddleSliceClusters,'stable'); % identify cluster indice
    
    % brainstem outline right
    cluster1 = positionLR_sorted(MiddleSliceClusters==ordersLR(4),:);
    transformed_landmark_coords{subj,1}.MidSlice.OutlineBrainstem_right = [median(cluster1(:,1)) median(cluster1(:,2)) slices(2)]; %ceil( median(cluster1) );
    
    % LC right
    cluster2 = positionLR_sorted(MiddleSliceClusters==ordersLR(3),:);
    transformed_landmark_coords{subj,1}.MidSlice.LC_right = [median(cluster2(:,1)) median(cluster2(:,2)) slices(2)];
    
    % LC left
    cluster3 = positionLR_sorted(MiddleSliceClusters==ordersLR(2),:);
    transformed_landmark_coords{subj,1}.MidSlice.LC_left = [median(cluster3(:,1)) median(cluster3(:,2)) slices(2)];
    
    % brainstem outline left
    cluster4 = positionLR_sorted(MiddleSliceClusters==ordersLR(1),:);
    transformed_landmark_coords{subj,1}.MidSlice.OutlineBrainstem_left = [median(cluster4(:,1)) median(cluster4(:,2)) slices(2)];
    
    
    %% the third slice, with one landmark: bottom window brainstem
    
    clear idx3 positionSlice positionLR_sorted cluster1 cluster2 cluster3 BottomSliceClusters
    
    transformed_landmark_coords{subj,1}.BottomSlice = [];
    
    positionSlice = tmp_sorted(find(tmp_sorted(:,3)==slices(3)),:);
    [~,idx3] = sort(positionSlice(:,1)); positionLR_sorted = positionSlice(idx3,:);
    
    % brainstem bottom
    transformed_landmark_coords{subj,1}.BottomSlice.BottomBrainstem = [median(positionLR_sorted(:,1)) median(positionLR_sorted(:,2)) slices(3)]; %ceil( median(positionLR_sorted) );
    
    fprintf('\n subject %s done\n',IDs{subj})
    
end;clear subj
 
disp('coordinates collected')
 
%% calculate distances in single subject level
 
varnames = {'NucRuber_left';'NucRuber_right';'TopBrainstem';'OutlineBrainstem_left';'OutlineBrainstem_right';...
    'LC_left';'LC_right';'BottomBrainstem'};
for v1=1:length(varnames)
   eval([varnames{v1} '=[];']) 
end
 
Distances_indv=[]; 
for subj=1:length(IDs)
    
    Distances_indv{subj,1}.NucRuber=[];
    Distances_indv{subj,1}.NucRuber.left = sum((transformed_landmark_coords{subj,1}.TopSlice.NucRuber_left-MNI.NucRuber.left).^2).^0.5;
    Distances_indv{subj,1}.NucRuber.right = sum((transformed_landmark_coords{subj,1}.TopSlice.NucRuber_right-MNI.NucRuber.right).^2).^0.5;
    
    Distances_indv{subj,1}.TopBrainstem = sum((transformed_landmark_coords{subj,1}.TopSlice.TopBrainstem-MNI.TopBrainstem).^2).^0.5;
    
    Distances_indv{subj,1}.OutlineBrainstem=[];
    Distances_indv{subj,1}.OutlineBrainstem.left=sum((transformed_landmark_coords{subj,1}.MidSlice.OutlineBrainstem_left-MNI.OutlineBrainstem.left).^2).^0.5;
    Distances_indv{subj,1}.OutlineBrainstem.right=sum((transformed_landmark_coords{subj,1}.MidSlice.OutlineBrainstem_right-MNI.OutlineBrainstem.right).^2).^0.5;
    
    Distances_indv{subj,1}.LC=[];
    Distances_indv{subj,1}.LC.left=sum((transformed_landmark_coords{subj,1}.MidSlice.LC_left-MNI.LC.left).^2).^0.5;
    Distances_indv{subj,1}.LC.right=sum((transformed_landmark_coords{subj,1}.MidSlice.LC_right-MNI.LC.right).^2).^0.5;
    
    Distances_indv{subj,1}.BottomBrainstem = sum((transformed_landmark_coords{subj,1}.BottomSlice.BottomBrainstem-MNI.BottomBrainstem).^2).^0.5;
    
    % write as a table
    NucRuber_left(subj,1)=Distances_indv{subj,1}.NucRuber.left;
    NucRuber_right(subj,1)=Distances_indv{subj,1}.NucRuber.right;
    TopBrainstem(subj,1)=Distances_indv{subj,1}.TopBrainstem;
    OutlineBrainstem_left(subj,1)=Distances_indv{subj,1}.OutlineBrainstem.left;
    OutlineBrainstem_right(subj,1)=Distances_indv{subj,1}.OutlineBrainstem.right;
    LC_left(subj,1)=Distances_indv{subj,1}.LC.left;
    LC_right(subj,1)=Distances_indv{subj,1}.LC.right;
    BottomBrainstem(subj,1)=Distances_indv{subj,1}.BottomBrainstem;
    
end
 
IDcolumn = cell2mat(cellfun(@str2num, IDs,'UniformOutput',false));
 
T=table(NucRuber_left,NucRuber_right,TopBrainstem,OutlineBrainstem_left,OutlineBrainstem_right,...
    LC_left,LC_right,BottomBrainstem);
Distance_export =[IDcolumn',NucRuber_left,NucRuber_right,TopBrainstem,OutlineBrainstem_left,OutlineBrainstem_right,...
    LC_left,LC_right,BottomBrainstem];
 
save(['/Users/alex/Documents/ALEX_landmarks/EuclideanDistance_landmarks_FallypridePilot.mat'],'Distance_export')
 
disp('distance calc done')

%% calculate interrator agreement

load('/Users/alex/Dropbox/paperwriting/coreg/data/ID_pilot.mat')

clear transformed_landmark_coords
ML=load('/Users/alex/Documents/ALEX_landmarks/LandmarkCoordinates_pilot_ML.mat');
AY=load('/Users/alex/Desktop/rt_epis_new/LandmarkCoordinates_pilot_AY.mat');

% calculate agreements per subject per landmark
IR=[];
for id=1:length(IDs)
    
    
    % --------- top slice --------- %
    
    % top slice
    clear ay ml
    ay=AY.transformed_landmark_coords{id,1}.TopSlice;
    ml=ML.transformed_landmark_coords{id,1}.TopSlice;
    IR{id,1}.NucRuber_left      =sum(ay.NucRuber_left(1:2)==ml.NucRuber_left(1:2));
    IR{id,1}.NucRuber_right     =sum(ay.NucRuber_right(1:2)==ml.NucRuber_right(1:2));
    IR{id,1}.TopBrainstem       =sum(ay.TopBrainstem(1:2)==ml.TopBrainstem(1:2));
    
    NucRuber_left(id,1)         =sum(ay.NucRuber_left(1:2)==ml.NucRuber_left(1:2));
    NucRuber_right(id,1)        =sum(ay.NucRuber_right(1:2)==ml.NucRuber_right(1:2));
    TopBrainstem(id,1)          =sum(ay.TopBrainstem(1:2)==ml.TopBrainstem(1:2));
    
    
    % middle slice
    clear ay ml
    ay=AY.transformed_landmark_coords{id,1}.MidSlice;
    ml=ML.transformed_landmark_coords{id,1}.MidSlice;
    IR{id,1}.OutlineBrainstem_left      =sum(ay.OutlineBrainstem_left(1:2)==ml.OutlineBrainstem_left(1:2));
    IR{id,1}.OutlineBrainstem_right     =sum(ay.OutlineBrainstem_right(1:2)==ml.OutlineBrainstem_right(1:2));
    IR{id,1}.LC_left                    =sum(ay.LC_left(1:2)==ml.LC_left(1:2));
    IR{id,1}.LC_right                   =sum(ay.LC_right(1:2)==ml.LC_right(1:2));
    
    OutlineBrainstem_left(id,1)       =sum(ay.OutlineBrainstem_left(1:2)==ml.OutlineBrainstem_left(1:2));
    OutlineBrainstem_right(id,1)      =sum(ay.OutlineBrainstem_right(1:2)==ml.OutlineBrainstem_right(1:2));
    LC_left(id,1)                     =sum(ay.LC_left(1:2)==ml.LC_left(1:2));
    LC_right(id,1)                    =sum(ay.LC_right(1:2)==ml.LC_right(1:2));
    
    % bottom slice
    clear ay ml
    ay=AY.transformed_landmark_coords{id,1}.BottomSlice;
    ml=ML.transformed_landmark_coords{id,1}.BottomSlice;
    IR{id,1}.BottomBrainstem        =sum(ay.BottomBrainstem(1:2)==ml.BottomBrainstem(1:2));
    BottomBrainstem(id,1)           =sum(ay.BottomBrainstem(1:2)==ml.BottomBrainstem(1:2));
    
end

%% calculate DICE score

DICE=[];

DICE.NucRuber_left=(sum(NucRuber_left)*2)/(2*2*numel(IDs));
DICE.NucRuber_right=(sum(NucRuber_right)*2)/(2*2*numel(IDs));
DICE.TopBrainstem=(sum(TopBrainstem)*2)/(2*2*numel(IDs));

DICE.OutlineBrainstem_left=(sum(OutlineBrainstem_left)*2)/(2*2*numel(IDs));
DICE.OutlineBrainstem_right=(sum(OutlineBrainstem_right)*2)/(2*2*numel(IDs));
DICE.LC_left=(sum(LC_left)*2)/(2*2*numel(IDs));
DICE.LC_right=(sum(LC_right)*2)/(2*2*numel(IDs));

DICE.BottomBrainstem=(sum(BottomBrainstem)*2)/(2*2*numel(IDs));
