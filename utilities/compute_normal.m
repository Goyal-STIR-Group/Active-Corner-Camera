function [ objects ] = compute_normal( objects, k )
nmz=@(x) x/norm(x); %Normalize vector


objects{k,5} = -nmz(cross(objects{k,3},objects{k,4})); % normal vector


if norm(objects{k,2}+objects{k,5})>norm(objects{k,2})
    objects{k,5} = -objects{k,5};
end
end

