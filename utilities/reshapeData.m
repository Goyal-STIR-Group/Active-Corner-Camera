function [dataMat] = reshapeData(dataVec,cam_pixel_dim,num_time_bins)

dataMat = zeros(cam_pixel_dim^2,num_time_bins);
pullIndex = 1;
for i = 1:cam_pixel_dim^2
%     for j = 1:cam_pixel_dim
        dataMat(i,:) = dataVec(pullIndex:pullIndex+num_time_bins-1);
        pullIndex = pullIndex+num_time_bins;
%     end
end
end

