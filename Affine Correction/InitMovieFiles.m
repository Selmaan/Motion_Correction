%% Change This: Create variable / filenames names and location
num_files = 10;
nameOffset = 0;
mouse_name = 'AM115';
session_name = '1_3x';
view_name = '119Sub';
slice_name = '001';
movie_file = sprintf('%s_%s_%s_%s',mouse_name,session_name,view_name,slice_name);
ICPC_file = sprintf('%s_ICPC',movie_file);
tiffPath = 'E:\Data\Ari\AM115';
filepath = 'E:\Data\Corrected Files\AM115\4 Segs';


for j=1:num_files
    correct_filenames{j} = [tiffPath '\' sprintf('%s_%s_%s_%s_%.3d_red.tif',...
        mouse_name,session_name,view_name,slice_name,j+nameOffset)];
    apply_filenames{j} = [tiffPath '\' sprintf('%s_%s_%s_%s_%.3d_green.tif',...
        mouse_name,session_name,view_name,slice_name,j+nameOffset)];
end

%% Do Not Change This: Initialize File Variables
cd(filepath),
MovFile = matfile([movie_file '.mat'],'Writable',true);
MovFile.movie_mask = [];
MovFile.acqFrames=[];
MovFile.cated_xShift = [];
MovFile.cated_yShift = [];
MovFile.acqRef = zeros(0,0,0,'single');
MovFile.correct_filenames = correct_filenames;
MovFile.apply_filenames = apply_filenames;