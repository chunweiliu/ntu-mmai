function cbir( q_set )
%% create feature vectors

%feature_vectors_color( );
%feature_vectors_texture( );
%feature_vectors_color_texture( );

%% build a tree and read the index file

fid_fv_path = fopen( 'database\features\fv_path.txt', 'r' );
fv_path = [];
while 1
    tline = fgetl( fid_fv_path );
    if ~ischar( tline )
        break;
    end
    c = cellstr(tline);         % make string to cell let array to save them
    fv_path = [ fv_path; c];    % cannot deal with the different size row
end
fclose( fid_fv_path );

query_color_texture( fv_path, q_set );

return;

%% read all *.jpg files in the folder path
function feature_vectors_color( )

path = 'dataset\';
folders = {'dress'};
type = '\*.jpg';

fv_path = 'database\features\fv_path.txt';
fid_fv_path = fopen( fv_path, 'w+' );

fvc_list = [];
for n = 1 : length( folders )

    filetypes = strcat(path, folders(n), type); % query all .jpg file
    link = cell2mat( filetypes );
    files = dir( link );

    for m = 1 : length( files )
        filepath = strcat(path, folders(n), '\', files(m).name); % after get the file name, we add the file path here
        link = cell2mat( filepath );
        img = imread( link );

        % process the image here
        fvc = fv_color( img );     % 3-d array
        fvc_list = [fvc_list; fvc];

        % save the fv path
        fprintf( fid_fv_path, '%s\n', link);
    end
end
fclose( fid_fv_path );
fvc_path = 'database\features\fvc_list.mat';
save( fvc_path, 'fvc_list' );
return;

%%
function feature_vectors_texture( )

path = 'dataset\';
folders = {'dress'};
type = '\*.jpg';

fvt_list = [];
for n = 1 : length( folders )

    filetypes = strcat(path, folders(n), type); % query all .jpg file
    link = cell2mat( filetypes );
    files = dir( link );

    for m = 1 : length( files )
        filepath = strcat(path, folders(n), '\', files(m).name); % after get the file name, we add the file path here
        link = cell2mat( filepath );
        img = imread( link );

        % process the image here
        fvt = fv_texture( img );
        fvt_list = [fvt_list; fvt];

    end
end
fvt_path = 'database\features\fvt_list.mat';
save( fvt_path, 'fvt_list' );
return;

%% average of color and texture similarity
function feature_vectors_color_texture(  )
fvc_list_struct = load( 'database\features\fvc_list.mat' );
fvc_list_cell = struct2cell( fvc_list_struct );
fvc_list = cell2mat( fvc_list_cell );

fvt_list_struct = load( 'database\features\fvt_list.mat' );
fvt_list_cell = struct2cell( fvt_list_struct );
fvt_list = cell2mat( fvt_list_cell );


[r, c] = size( fvt_list );
fvm_list = [];
for n = 1:c
    if mod(n, 2) == 1
        m = ceil(n/2);
        v = fvt_list(:, n) + 10.* fvc_list(:, m);
    else
        v = fvt_list(:, n);
    end
    fvm_list = [fvm_list, v];
end
fvm_path = 'database\features\fvm_list.mat';
save( fvm_path, 'fvm_list' );
return;

%%
function query_color_texture( fv_path, q_set )
fvm_list_struct = load( 'database\features\fvm_list.mat' );
fvm_list_cell = struct2cell( fvm_list_struct );
fvm_list = cell2mat( fvm_list_cell );
[row, col] = size( fvm_list );
tree_fvt = kdtree( fvm_list(:, 1:col) );

img = imread( q_set );
cell_path = im_query( tree_fvt, fv_path, img, @fv_ct );
show_image( cell_path );
return;

%% produce the color feature vector
function fv = fv_color( img )
img3 = if_gray( img );
img_hsv = rgb2hsv( img3 );
bins_list = [16, 4, 4];

fv = [];
for k = 1 : 3
    img_tmp = img_hsv( :, :, k );
    bins = bins_list( k );
    [n, b] = imhist( img_tmp, bins );

    [r, c] = size( img_tmp );
    norm_n = n ./ (r*c);
    fv = [fv, norm_n'];          % row vector
end
return;

%% gabor filter make fv_texture
function fv = fv_texture( img )
nscale          = 4;        %Number of wavelet scales.
norient         = 6;        %Number of filter orientations.
minWaveLength   = 3;        %Wavelength of smallest scale filter.
mult            = 2;        %Scaling factor between successive filters.
sigmaOnf        = 0.65;     %Ratio of the standard deviation of the
                            %Gaussian describing the log Gabor filter's
                            %transfer function in the frequency domain
                            %to the filter center frequency.
dThetaOnSigma   = 1.5;      %Ratio of angular interval between filter
                            %orientations and the standard deviation of
                            %the angular Gaussian function used to
                            %construct filters in the freq. plane.

img3 = if_gray( img );
img_gray = rgb2gray( img3 );
leng = length( img_gray );
img_re = imresize( img_gray, 500/leng );

gacell = gaborconvolve( img_re, nscale, norient, minWaveLength, ...
    mult, sigmaOnf, dThetaOnSigma);
fv = [];
for n = 1: nscale
    for m = 1: norient
        gamat = cell2mat( gacell(n, m) );
        gamean = mean( mean( abs( gamat ) ) );
        gastd = std( std( abs( gamat ) ) );
        fvone = [gamean, gastd];  % row vector 1 * 48
        fv = [fv, fvone];
    end
end
return;

%%
function fv = fv_ct( img )
fvc = fv_color( img );
fvt = fv_texture( img );

[r, c] = size( fvt );
fv = [];
for n = 1:c
    if mod(n, 2) == 1
        m = ceil(n/2);
        v = fvt(1, n) + 10*fvc(1, m);
    else
        v = fvt(1, n);
    end
    fv = [fv, v];
end
return;

%% gray level
function re_img = if_gray( img )
[a, b, c] = size( img );
if ( c == 1 )
    re_img = [];
    for n = 1 : 3
        re_img(:, :, n) = img(:, :, 1);
    end
else
    re_img = img;
end
return;

%% im_query function
function cell_path = im_query( tree, path, img, fv )
q = fv( img );

indx_color = kdtree_closestpoint(tree, q);

cell_path = path( indx_color, 1 );
return;

%% show image function
function show_image( cell_path )
fo_path = 'outputs\query_result.txt';
fid = fopen( fo_path, 'w' );
%show_path = cell2mat( cell_path );
A=(char(regexp(cell_path,'\d+','match','once')));
fwrite(fid, A);
%s = sprintf('%s', show_path);
%fwrite(fid, s);
fclose(fid);
return;
