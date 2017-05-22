%This will generate images and video in subfolder called 'images'!
function sign = geospl(filename) % Designed for data downloaded from NCEP, where data is lon*lat*level.
lat = ncread(filename, 'lat');
lon = ncread(filename, 'lon');
slp = ncread(filename, 'slp');
latlim = double([min(lat) max(lat)]); % Range of latitude and longitude.
lonlim = double([min(lon) max(lon)]);
dim = size(slp);
num_obs = dim(3); % Number of observations.
imageNames = cell(num_obs,1);
for i = 1:num_obs
    exp = transpose(slp(:,:,i)); % Need transpose for lat*lon.
    R = georasterref('RasterSize', size(exp), 'Latlim', latlim, 'Lonlim', lonlim) % Set up raster reference from dataset.
    figure('Visible','off','Color','k'); 
    worldmap(latlim, lonlim); % Set up world map.
    load geoid
    geoshow(geoid, geoidrefvec, 'DisplayType', 'texturemap');
    geoshow('landareas.shp', 'FaceColor', [0.15 0.5 0.15]); % Add land areas and coastlines.
    geoshow(exp, R, 'DisplayType', 'texturemap');
    mesh(peaks);
    caxis([min(min(exp)) max(max(exp))])
    colormap(parula);
    figure_nm = strcat(pwd,'/images/GeoSlp_', num2str(i),'.bmp');
    imageNames(i) = cellstr(figure_nm); % Save image directory.
    saveas(gcf,figure_nm); % Save images.
end;
outputVideo = VideoWriter(strcat(pwd,'/SeaPressure.avi'));
outputVideo.FrameRate = 5; % Set FPS manually.
open(outputVideo);
for i = 1:num_obs
   img = imread(char(imageNames(i))); % Add images sequentially.
   writeVideo(outputVideo,img);
end
close(outputVideo);
sign = 'Done';
    

