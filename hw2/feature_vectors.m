function feature_vectors( )
%% read all *.jpg files in the folder path
    path = 'dataset\';
    % use {} to concat the string list, [] to concat the chars
    %folders = {'bag', 'comic', 'dog', 'flowers', 'food', 'instrument', ...
    %           'money','painting', 'starbucks', 'taipei101', 'tool', 'tracffic'};
    folders = {'bag'};
    type = '\*.jpg';
    for n = 1 : length( folders )
        
        filetypes = strcat(path, folders(n), type); % query all .jpg file
        link = cell2mat( filetypes );
        files = dir( link );
        
        for m = 1 : length( files )
            filepath = strcat(path, folders(n), '\', files(m).name); % after get the file name, we add the file path here
            link = cell2mat( filepath );
            img = imread( link );
            % process the image here
            img3 = if_gray( img );
            fv = fv_color( img3 );
            
            % save the fv
            savefile = strcat(path, folders(n), '\', 'fv_color_',files(m).name); % .jpg should not in your file name !!!!
            link = cell2mat( savefile );
            fid = fopen( link, 'w' );
            fwrite( fid, '%f', fv);
            fclose( fid );
        end
        
    end
return;

%% produce the color feature vector
function fv = fv_color( img )
    hsv_img = rgb2hsv( img );
    h = hsv_img( :, :, 1 );         % hue color name
    s = hsv_img( :, :, 2 );         % saturation [gray, brightness]
    v = hsv_img( :, :, 3 );         % value [dark, light]
    hist_list = [h, s, v];
    bins_list = [16, 4, 4];
    
    fv = [];
    for k = 1 : 3
        img = hist_list( k );
        bins = bins_list( k );
        [n, b] = hist( img, bins );
        
        [r, c] = size( img );
        norm_n = n ./ (r*c/bins);
        fv = [fv; norm_n'];
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