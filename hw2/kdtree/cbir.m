function cbir( flag )
%% create feature vectors
    if( flag == 0 )
        feature_vectors( );

%% build a tree and read the index file
    elseif( flag == 1 )   
        % mex kdtree_build.cpp
        mex kdtree_k_nearest_neighbors.cpp
        % disp('compiled.');
        
        fid_fv_path = fopen( 'output\fv_path.txt', 'r' );
        fv_path = [];
        while 1
            tline = fgetl( fid_fv_path );
            if ~ischar( tline )
                break;
            end
            c = cellstr(tline);         % make string to cell let array to save them
            fv_path = [ fv_path; c];  % cannot deal with the different size row
        end 
        fclose( fid_fv_path );

        fvc_list_struct = load( 'output\fvc_list.mat' );
        fvc_list_cell = struct2cell( fvc_list_struct );
        fvc_list = cell2mat( fvc_list_cell );
        tree_fvc = kdtree_build( fvc_list(:, 1:24) );        
        
        fvt_list_struct = load( 'output\fvt_list.mat' );
        fvt_list_cell = struct2cell( fvt_list_struct );
        fvt_list = cell2mat( fvt_list_cell );
        tree_fvt = kdtree_build( fvt_list(:, 1:48) );
%% strat query and save the result
        q_times = 0;
        q_set = {'dataset\00queries00\money\q_money01.jpg', ...
            'dataset\00queries00\starbucks\q_starbucks02.jpg', ...
            'dataset\00queries00\sunflower\q_sunflower02.jpg', ...
            'dataset\00queries00\taipei101\q_taipei03.jpg'};
        for cq_path = q_set
        %while 1
            q_times = q_times + 1;
            %q_path = input('Please input the file path of image:\n', 's' );
            %if isempty( query_path )
            %    break;
            %end
            q_path = cell2mat( cq_path );
            img = imread( q_path );
            cell_path_c = im_query( tree_fvc, fv_path, img, 10 );
            show_image( q_times, cell_path_c, 1 );
            
            cell_path_t = im_query( tree_fvt, fv_path, img, 10 );
            show_image( q_times, cell_path_t, 2 );
        end
    end
return;

%% read all *.jpg files in the folder path
function feature_vectors( )

    path = 'dataset\';
    % use {} to concat the string list, [] to concat the chars
    %folders = {'bag', 'comic', 'dog', 'flowers', 'food', 'instrument', ...
    %    'money', 'painting', 'starbucks', 'taipei101', 'tool', 'tracffic'};
    folders = {'bag', 'comic'};
    type = '\*.jpg';
    
    fv_path = 'output\fv_path.txt';
    fid_fv_path = fopen( fv_path, 'w+' );
    
    fvc_list = [];
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
            fvc = fv_color( img );     % 3-d array
            fvc_list = [fvc_list; fvc];
            
            fvt = fv_texture( img );
            fvt_list = [fvt_list; fvt];
            
            % save the fv path           
            fprintf( fid_fv_path, '%s\n', link);
        end
    end
    fclose( fid_fv_path );
    savepathc = 'output\fvc_list.mat';
    save( savepathc, 'fvc_list' );
    savepatht = 'output\fvt_list.mat';
    save( svaepatht, 'fvt_list' );
return;

%% produce the color feature vector
function fv = fv_color( img )
    img3 = if_gray( img );
    img_hsv = rgb2hsv( img3 );
    %h = img_hsv( :, :, 1 );         % hue color name
    %s = img_hsv( :, :, 2 );         % saturation [gray, brightness]
    %v = img_hsv( :, :, 3 );         % value [dark, light]
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
    nscale          = 4;      %Number of wavelet scales.
    norient         = 6;      %Number of filter orientations.
    minWaveLength   = 3;      %Wavelength of smallest scale filter.
    mult            = 2;      %Scaling factor between successive filters.
    sigmaOnf        = 0.65;   %Ratio of the standard deviation of the
                              %Gaussian describing the log Gabor filter's
                              %transfer function in the frequency domain
                              %to the filter center frequency. 
    dThetaOnSigma   = 1.5;    %Ratio of angular interval between filter
                              %orientations and the standard deviation of
                              %the angular Gaussian function used to
                              %construct filters in the freq. plane.
    img3 = if_gray( img );
    img_gray = rgb2gray( img3 );
    gacell = gaborconvolve( img_gray, nscale, norient, minWaveLength, ...
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
function cell_path = im_query( tree, path, img, k )    
    q = fv_color( img );
    q = q';

    indx_color = kdtree_k_nearest_neighbors(tree, q, k);      

    cell_path = path( indx_color, 1 );    
return;

%% show image function
function show_image( q_times, cell_path, k )
    for n = 1:length( cell_path )
        
        show_path = cell2mat( cell_path(n,1) );
        
        save_time = num2str( q_times );
        save_num = num2str( n );
        if k == 1
            save_path = strcat('query_results_color\', save_time, '_query_', save_num, '.png');
        elseif k == 2
            save_path = strcat('query_results_texture\', save_time, '_query_', save_num, '.png');
        end
        system( ['copy ' show_path ' ' save_path] );
    end
return;

%% gabor filter
% function gb=gabor_fn(sigma_x, sigma_y, theta, lambda, psi, gamma)
%  
%     sz_x = fix(6*sigma_x);
%     if mod(sz_x,2)==0, sz_x=sz_x+1;end
% 
%     sz_y = fix(6*sigma_y);
%     if mod(sz_y,2)==0, sz_y=sz_y+1;end
% 
%     [x y]=meshgrid(-fix(sz_x/2):fix(sz_x/2),fix(-sz_y/2):fix(sz_y/2));
% 
%     % Rotation 
%     x_theta=x*cos(theta)+y*sin(theta);
%     y_theta=-x*sin(theta)+y*cos(theta);
% 
%     gb=exp(-.5*(x_theta.^2/sigma_x^2+gamma^2*y_theta.^2/sigma_y^2)).*cos(2*pi/lambda*x_theta+psi);
% return;

