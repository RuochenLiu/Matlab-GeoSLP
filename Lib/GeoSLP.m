function sign = GeoSLP(nc_name, FPS)
%
% GeoSLP is designed to transform and plot sea level pressure data on
% world map to create a video. GeoSLP reads information of latitude, longtitude,
% sea level pressure matrix from NC files downloaded from NCEP website, then 
% converting them into mat files. Use ranges of latitude and longitude to create a
% part of world map with coastlines for plotting. Plotting is based on raster
% reference data and geoshow function.
%
% INPUT 
%   'nc_name' is the filename of sea level pressure NC file, which is downloaded
%   from NCEP website.
%
%   'FPS' is the value of fps of video to generate, with default 5.
% 		
% OUTPUT 
%   SeaLevelPressure_Info.mat, which stores information of latitude, longitude, 
%   sea level pressure values, will be generated under the working directory.
%
%   SeaPressure.avi, which stores geo-plot images, will be generated under 
%   the working directory.
% 
% Examples:
%   GeoSLP('filename.nc', 10);
%
%   GeoSLP('filename.nc') under default fps;
%
%
%% AUTHOR    : Ruochen Liu
%% DATE     : 17-May-2017
%% Revision : 1.00
%% DEVELOPED : R2016a
%% FILENAME  : GeoSLP.m
%
    if  nargin < 2
        FPS = 5; % Defalut fps value.
    else
        FPS = FPS;
    end

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
    
    t1 = datetime(2012,4,6);
	t2 = datetime(2014,12,31);
	dateSeq = t1:t2;
    
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
        title(datestr(dateSeq(i)), 'FontSize', 12);

		fig = gcf;
		fig.Color = 'white';
		F = getframe(gcf);
		writeVideo(outputVideo, F);
        clf;
	end

	close(outputVideo);
	disp('Complete');

end
