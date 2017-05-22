function sign = GeoSPL(nc_name, FPS)
%
% GeoSPL is designed to transform and plot sea pressure level data on
% world map to create a video. GeoSPL reads information of latitude, longtitude,
% sea pressure level matrix from NC files downloaded from NCEP website, then 
% converting them into mat files. Use ranges of latitude and longitude to create a
% part of world map with coastlines for plotting. Plotting is based on raster
% reference data and geoshow function.
%
% INPUT 
%   'nc_name' is the filename of sea pressure level NC file, which is downloaded
%   from NCEP website.
%
%   'FPS' is the value of fps of video to generate, with default 5.
% 		
% OUTPUT 
%   SeaPressureLevel_Info.mat, which stores information of latitude, longitude, 
%   sea pressure level values, will be generated under the working directory.
%
%   SeaPressure.avi, which stores geo-plot images, will be generated under 
%   the working directory.
% 
% Examples:
%   GeoSPL('filename.nc', 10);
%
%   GeoSPL('filename.nc') under default fps;
%
%
%% AUTHOR    : Ruochen Liu
%% DATE     : 17-May-2017
%% Revision : 1.00
%% DEVELOPED : R2016a
%% FILENAME  : GeoSPL.m
%
    if  nargin < 2,
        FPS = 5; % Defalut fps value.

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
    mat_name = 'SeaPressureLevel_Info.mat';
    save(mat_name, 'Dimension', 'Latitude', 'Longitude', 'SeaPressure'); % Save data as mat file.
    
    
    
	%imageNames = cell(num_obs,1);
	%imDirName = strcat(pwd,'/images/');
	%if exist(imDirName) ~= 7,  
	%	mkdir( imDirName );  
	%end

	crange = prctile(SeaPressure(:) , [1 99]);

	outputVideo = VideoWriter(strcat(pwd,'/SeaPressure.avi'));
	outputVideo.FrameRate = FPS; % Set FPS.
	open(outputVideo);

	for i = 1:num_obs,

		disp( strcat( 'Generating frame # ', num2str(i) ) );

		exp = SeaPressure(:,:,i);
		% FIXME Check the following to see if it does what I assumed it does.
		R = georasterref('RasterSize', size(exp), 'Latlim', latlim, 'Lonlim', lonlim); % Set up raster reference from dataset.
		figure('Visible','off','Color','k'); 
		worldmap(latlim, lonlim); % Set up world map.
		load geoid
		geoshow(geoid, geoidrefvec, 'DisplayType', 'texturemap');
		geoshow('landareas.shp', 'FaceColor', [0.15 0.5 0.15]); % Add land areas and coastlines.
		geoshow(exp, R, 'DisplayType', 'texturemap');
		mesh(peaks);
		caxis(crange);
		colormap(parula);

		fig = gcf;
		fig.Color = 'white';
		F = getframe(gcf);
		writeVideo(outputVideo, F);
	end

	close(outputVideo);
	sign = 'Done';

end
