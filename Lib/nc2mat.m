function nc2mat(nc_name)
lat = ncread(nc_name, 'lat');
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
end


