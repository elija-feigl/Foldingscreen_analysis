function imageData = load_gel_image(varargin)
%% Loads gel image data
%   INPUTS: 
%   'data_dir' (optional parameter) = initial directory, where data is stored
%   'n_images' (optional parameter) = number of images to load
%   'gel_format' (optional parameter) = pass 'on' to load .gel files
%   OUTPUT:
%   imageData struct with .images .pathnames .filenames .nrImages
%   .nrImages is number of images
%   .images is cell array {nr_image} of image data
%   .pathnames is cell array {nr_image} of pathnames to each image
%   .filenames is cell array {nr_image} of filenames of each image
%
% Example: myImageData = load_gel_image('data_dir', '~/Documents/Images');

%% parse input
p = inputParser;
default_data_dir = userpath; % default data directory is userpath
default_data_dir=default_data_dir(1:end-1);
dafault_paths_to_images = {}; 
% optional paramters

addParameter(p,'data_dir',default_data_dir, @isstr);
% n_images, default is -1
addParameter(p,'n_images', -1 ,@isnumeric);
% load files with .gel ending
addParameter(p,'gel_format', 'off', @(x) strcmp(x, 'on'));

addParameter(p,'paths_to_images',dafault_paths_to_images, @iscell);


parse(p,  varargin{:});
% default data location
data_dir = p.Results.data_dir;
% set number of images
nrImages = p.Results.n_images;
% .gel data format loading on/off
gel_format_bool = strcmp(p.Results.gel_format, 'on');
% path to images
paths_to_images = p.Results.paths_to_images;

%% select image data
if isempty(paths_to_images)
    init_path = cd; %remember initial/current path

    if nrImages <= 0 % ask how many images user wants to load
        temp = inputdlg({'How many images (channels) do you want to load:'}, 'How many images (channels) do you want to load?', 1, {'1'});
        nrImages = str2double(temp(1));
    end

    filenames = cell(nrImages, 1);
    pathnames = cell(nrImages, 1);

    lastDirectory = data_dir;
    for i=1:nrImages
        cd(lastDirectory)
        if gel_format_bool
            [filenames{i}, pathnames{i}]=uigetfile('*.gel','Select image:');
        else
            [filenames{i}, pathnames{i}]=uigetfile('*.tif','Select image:');
        end
        lastDirectory = pathnames{i};
    end
    cd(init_path) % cd to initial directory
else
    nrImages = length(length(paths_to_images));
    filenames = cell(nrImages, 1);
    pathnames = cell(nrImages, 1);
    %compute pathnames and 
    for i=1:nrImages
        [filepath,name,ext] = fileparts(paths_to_images{i});
        pathnames{i} = filepath;
        filenames{i} = [name ext];
    end
end
%% load image data
images = cell(nrImages, 1);

for i=1:nrImages
    images{i} = double(imread([pathnames{i} filesep filenames{i}]));             %load image data  
end

%% create imageData structure, return imageData structure

imageData=struct('images',{images},'pathnames',{pathnames},'filenames',{filenames},'nrImages',nrImages);