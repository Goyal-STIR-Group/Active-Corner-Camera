function [vecXYZ] = rsu2xyz(R,o,vecRSU,npoint)

% rotate back
vecRSU_rotateBack = R'*vecRSU';
vecXYZ = vecRSU_rotateBack' + repmat(o,npoint,1);

end

