%%  Visual check for transformed images
%   make a movie of transformed EPIs/T1s to see whether there's any obvious
%   misalignment
 
%% preparation
 
clear;clc
 
% set paths
path_root = '/path/to/your/transformed/images/'; % where are the files?
path_save = '/path/where/the/video/will/be/saved/'; % where is the video going to be saved?
 
% set variables
ids = [];
 
input_prompt = {'duration per frame in seconds (default=0)'; 'name of the file'; 'view(sagittal, coronal, axial)'; 'which slice should the video capture?'};
defaults     = {'0','null','sagittal','96'};
input_answer = inputdlg(input_prompt, 'specify the properties of the video', 1, defaults);
 
vidframerate = str2num(input_answer{1,1});
vidformat    = '.avi';
vidname      = [input_answer{2,1} vidformat];
vidview      = input_answer{3,1};
slicenum     = str2num(input_answer{4,1});
 
cd(path_save)
v            = VideoWriter(vidname);
if vidframerate ~= 0
    v.FrameRate  = 1/vidframerate;
elseif vidframerate == 0
    v.FrameRate  = 30; % it's super fast
end
 
%% start recording
 
open(v) % open the file
 
cc = 0; figure % open a sketchbook
for i = 1:length(ids)        
        name = num2str(ids(i)); disp(num2str(ids(i)))
        data = spm_read_vols(spm_vol([path_root name '/data/transformed_image.nii']));

        % position the slice - dimensions are: (x=sagittal, y=coronal, z=axial)
        if strcmpi(vidview,'sagittal')
            dum = squeeze(data(slicenum,:,:)); % sagittal
        elseif strcmpi(vidview,'coronal')
            dum = squeeze(data(:,slicenum,:)); % coronal view
        elseif strcmpi(vidview,'axial')
            dum = squeeze(data(:,:,slicenum)); % axial
        else
            error('check your view input')
        end
        
        dum = rot90(dum);% t1WB on template
%     dum = zscore(dum);
 
        imagesc(dum);
        title(num2str(ids(i)))
        cc = cc+1;
        % Store the frame
        M(cc)=getframe(gcf); % leaving gcf out crops the frame in the movie.
             
end
 
writeVideo(v,M) % export the video
close(v)
close all;
