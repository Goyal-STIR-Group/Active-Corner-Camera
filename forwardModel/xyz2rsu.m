function [vec_rsu] = xyz2rsu(R,o,vec_xyz,npoint)
    vec_xyz_translate = vec_xyz - repmat(o,npoint,1);
    vec_rsu = R*vec_xyz_translate';
    vec_rsu = vec_rsu';
end

