%% LANDMARK PLOTTING: create 3D plots of created heatmaps, and compare
 
clc;clear;close all
 
% set env
path_transformed = '/path/to/your/individual/landmarks/';
load('/path/to/your/subject/IDs/IDs.mat') % load IDs
 
% load individual landmark images
transformed_landmarks_single=[];
for subj=1:length(IDs)
    transformed_landmarks_single{subj,1}=spm_read_vols(spm_vol([ path_transformed 'landmark_' IDs{subj} '.nii']));
end; clear subj
 
% load the landmark image drawn on MNI
MNIlandmark = spm_read_vols(spm_vol(['/path/to/the/MNI/landmark/MNI_landmark.nii'])); % change here accordingly to your file name
MNIlandmark(find(MNIlandmark<1))=NaN;
 
% write down the coordinates of MNI landmarks
MNI.NucRuber=[];
MNI.NucRuber.left=[102 113 71];
MNI.NucRuber.right=[92 113 71];
 
MNI.TopBrainstem=[97 102 71];
 
MNI.OutlineBrainstem=[];
MNI.OutlineBrainstem.left=[109 103 60];
MNI.OutlineBrainstem.right=[85 103 60];
 
MNI.LC=[];
MNI.LC.left=[100 97 60];
MNI.LC.right=[94 97 60];
 
MNI.BottomBrainstem=[97 94 50];
 
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
    clear idx3 positionSlice positionX_sorted cluster1 cluster2 cluster3
    
    positionSlice = tmp_sorted(find(tmp_sorted(:,3)==slices(1)),:); % find the points in the first slice
    [~,idx3] = sort(positionSlice(:,1)); positionX_sorted = positionSlice(idx3,:);
    
    % ruber right
    cluster1 = min(positionX_sorted(:,1));
    cluster_med1 = positionX_sorted(find(positionX_sorted(:,1)==cluster1+1),:);
    cluster_medpt{subj,1} = cluster_med1( find( cluster_med1(:,2)==ceil( median(cluster_med1(:,2))) ),: );
    
    % top curve
    cluster2 = min(positionX_sorted(positionX_sorted(:,1)>cluster1+2,1));
    cluster_med2 = positionX_sorted(find(positionX_sorted(:,1)==cluster2+1),:);
    cluster_medpt{subj,2} = cluster_med2( find( cluster_med2(:,2)==ceil( median(cluster_med2(:,2))) ),: );
    
    % ruber left
    cluster3 = min(positionX_sorted(positionX_sorted(:,1)>cluster2+2,1));
    cluster_med3 = positionX_sorted(find(positionX_sorted(:,1)==cluster3+1),:);
    cluster_medpt{subj,3} = cluster_med3( find( cluster_med3(:,2)==ceil( median(cluster_med3(:,2))) ),: );
    
    transformed_landmark_coords{subj,1}.TopSlice = [];
    transformed_landmark_coords{subj,1}.TopSlice.NucRuber_left = cluster_medpt{subj,3};
    transformed_landmark_coords{subj,1}.TopSlice.NucRuber_right = cluster_medpt{subj,1};
    transformed_landmark_coords{subj,1}.TopSlice.TopBrainstem = cluster_medpt{subj,2};
    clear cluster_medpt
    
    
    %% the second slice, with four landmarks: brainstem outline and LC
    
    clear idx3 positionY positionX_sorted cluster1 cluster2 cluster3
    
    positionSlice = tmp_sorted(find(tmp_sorted(:,3)==slices(2)),:);
    [~,idx3] = sort(positionSlice(:,1)); positionX_sorted = positionSlice(idx3,:);
    
    % brainstem outline right
    cluster1 = min(positionX_sorted(:,1));
    cluster_med1 = positionX_sorted(find(positionX_sorted(:,1)==cluster1+1),:);
    cluster_medpt{subj,1} = cluster_med1( find( cluster_med1(:,2)==ceil( median(cluster_med1(:,2))) ),: );
    
    % LC right
    cluster2 = min(positionX_sorted(positionX_sorted(:,1)>(cluster1+2),1));
    cluster_med2 = positionX_sorted(find(positionX_sorted(:,1)==cluster2+1),:);
    cluster_medpt{subj,2} = cluster_med2( find( cluster_med2(:,2)==ceil( median(cluster_med2(:,2))) ),: );
    
    % LC left
    cluster3 = min(positionX_sorted(positionX_sorted(:,1)>cluster2+2,1));
    cluster_med3 = positionX_sorted(find(positionX_sorted(:,1)==cluster3+1),:);
    cluster_medpt{subj,3} = cluster_med3( find( cluster_med3(:,2)==ceil( median(cluster_med3(:,2))) ),: );
    
    % brainstem outline left
    cluster4 = min(positionX_sorted(positionX_sorted(:,1)>cluster3+2,1));
    cluster_med4 = positionX_sorted(find(positionX_sorted(:,1)==cluster4+1),:);
    cluster_medpt{subj,4} = cluster_med4( find( cluster_med4(:,2)==ceil( median(cluster_med4(:,2))) ),: );
    
    transformed_landmark_coords{subj,1}.MidSlice = [];
    transformed_landmark_coords{subj,1}.MidSlice.OutlineBrainstem_left = cluster_medpt{subj,4};
    transformed_landmark_coords{subj,1}.MidSlice.OutlineBrainstem_right = cluster_medpt{subj,1};
    transformed_landmark_coords{subj,1}.MidSlice.LC_left = cluster_medpt{subj,3};
    transformed_landmark_coords{subj,1}.MidSlice.LC_right = cluster_medpt{subj,2};
    clear cluster_medpt
    
    
    %% the third slice, with one landmark: bottom window brainstem
    
    clear idx3 positionY positionX_sorted cluster1 cluster2 cluster3 cluster4
    
    positionSlice = tmp_sorted(find(tmp_sorted(:,3)==slices(3)),:);
    [~,idx3] = sort(positionSlice(:,1)); positionX_sorted = positionSlice(idx3,:);
    
    % brainstem bottom
    cluster1 = min(positionX_sorted(:,1));
    cluster_med1 = positionX_sorted(find(positionX_sorted(:,1)==cluster1+1),:);
    cluster_medpt{subj,1} = cluster_med1( find( cluster_med1(:,2)==ceil( median(cluster_med1(:,2))) ),: );
    
    transformed_landmark_coords{subj,1}.BottomSlice = [];
    transformed_landmark_coords{subj,1}.BottomSlice.BottomBrainstem = cluster_medpt{subj,1}
    clear cluster_medpt
    
    
    fprintf('\n subject %s done\n',IDs{subj})
    
end;clear subj
 
disp('coordinates collected')
 
%% calculate euclidian distances in single subject level
 
varnames = {'NucRuber_left';'NucRuber_right';'TopBrainstem';'OutlineBrainstem_left';'OutlineBrainstem_right';...
    'LC_left';'LC_right';'BottomBrainstem'};
for v1=1:length(varnames)
   eval([varnames{v1} '=[];']) 
end
 
EuclideanDistances_indv=[]; 
for subj=1:length(IDs)
    
    EuclideanDistances_indv{subj,1}.NucRuber=[];
    EuclideanDistances_indv{subj,1}.NucRuber.left = sum((transformed_landmark_coords{subj,1}.TopSlice.NucRuber_left-MNI.NucRuber.left).^2).^0.5;
    EuclideanDistances_indv{subj,1}.NucRuber.right = sum((transformed_landmark_coords{subj,1}.TopSlice.NucRuber_right-MNI.NucRuber.right).^2).^0.5;
    
    EuclideanDistances_indv{subj,1}.TopBrainstem = sum((transformed_landmark_coords{subj,1}.TopSlice.TopBrainstem-MNI.TopBrainstem).^2).^0.5;
    
    EuclideanDistances_indv{subj,1}.OutlineBrainstem=[];
    EuclideanDistances_indv{subj,1}.OutlineBrainstem.left=sum((transformed_landmark_coords{subj,1}.MidSlice.OutlineBrainstem_left-MNI.OutlineBrainstem.left).^2).^0.5;
    EuclideanDistances_indv{subj,1}.OutlineBrainstem.right=sum((transformed_landmark_coords{subj,1}.MidSlice.OutlineBrainstem_right-MNI.OutlineBrainstem.right).^2).^0.5;
    
    EuclideanDistances_indv{subj,1}.LC=[];
    EuclideanDistances_indv{subj,1}.LC.left=sum((transformed_landmark_coords{subj,1}.MidSlice.LC_left-MNI.LC.left).^2).^0.5;
    EuclideanDistances_indv{subj,1}.LC.right=sum((transformed_landmark_coords{subj,1}.MidSlice.LC_right-MNI.LC.right).^2).^0.5;
    
    EuclideanDistances_indv{subj,1}.BottomBrainstem = sum((transformed_landmark_coords{subj,1}.BottomSlice.BottomBrainstem-MNI.BottomBrainstem).^2).^0.5;
    
    % write as a table
    NucRuber_left(subj,1)=EuclideanDistances_indv{subj,1}.NucRuber.left;
    NucRuber_right(subj,1)=EuclideanDistances_indv{subj,1}.NucRuber.right;
    TopBrainstem(subj,1)=EuclideanDistances_indv{subj,1}.TopBrainstem;
    OutlineBrainstem_left(subj,1)=EuclideanDistances_indv{subj,1}.OutlineBrainstem.left;
    OutlineBrainstem_right(subj,1)=EuclideanDistances_indv{subj,1}.OutlineBrainstem.right;
    LC_left(subj,1)=EuclideanDistances_indv{subj,1}.LC.left;
    LC_right(subj,1)=EuclideanDistances_indv{subj,1}.LC.right;
    BottomBrainstem(subj,1)=EuclideanDistances_indv{subj,1}.BottomBrainstem;
    
end
 
spssID = cell2mat(cellfun(@str2num, IDs,'UniformOutput',false));
 
T=table(NucRuber_left,NucRuber_right,TopBrainstem,OutlineBrainstem_left,OutlineBrainstem_right,...
    LC_left,LC_right,BottomBrainstem);
ED_export =[spssID',NucRuber_left,NucRuber_right,TopBrainstem,OutlineBrainstem_left,OutlineBrainstem_right,...
    LC_left,LC_right,BottomBrainstem];
 
save(['/path/where/the/folder/the/data/will/be/saved/EuclideanDistance_landmarks.mat'],'ED_export')
 
disp('ED calc done')
