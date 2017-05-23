%%
nc_name = 'pressure.nc';

geospl(nc_name, 5);

lat = ncread(nc_name, 'lat'); % Read latitdue, longitude and pressure matrix from NC file.
lon = ncread(nc_name, 'lon');
SeaPressure = ncread(nc_name, 'slp');
SeaPressure = permute(SeaPressure, [2 1 3]);

latlim = double([min(lat) max(lat)]); % Range of latitude and longitude.
lonlim = double([min(lon) max(lon)]);
Dimension = size(SeaPressure); % Dimension of dataset latitude*longitude*pages.
num_obs = Dimension(3); % Number of observations.
Latitude = repmat(lat, 1, Dimension(2));
Longitude = repmat(transpose(lon), Dimension(1), 1);
mat_name = 'SeaLevelPressure_Info.mat';
save(mat_name, 'Dimension', 'Latitude', 'Longitude', 'SeaPressure'); % Save data as mat file.
%presseure_map = pcolor(lon, lat, exp);
%presseure_map.EdgeAlpha = 0;


%%
R = georasterref('RasterSize', size(exp), 'Latlim', [30 60], 'Lonlim', [150 230])
figure;
po = worldmap([30,60], [150,230]);
load geoid
geoshow(geoid, geoidrefvec, 'DisplayType', 'texturemap');
geoshow('landareas.shp', 'FaceColor', [0.15 0.5 0.15]);
geoshow(exp, R, 'DisplayType', 'texturemap');
mesh(peaks);
caxis([min(min(exp)) max(max(exp))])
colormap(parula);
colorbar;

%%
outputVideo = VideoWriter(strcat(pwd,'/SeaPressure.avi'));
outputVideo.FrameRate = 5;
open(outputVideo);
for i = 1:1000
    img = imread(char(imageNames(i)));
    writeVideo(outputVideo,img);
end
close(outputVideo);

%%
imageNames = cell(1000,1);
for i = 1:1000
    imageNames(i)=cellstr(strcat(pwd,'/images/GeoSlp_', num2str(i),'.bmp'));
end;