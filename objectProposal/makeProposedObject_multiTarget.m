function [proposedObjects] = makeProposedObject_multiTarget(phiBounds, heightVec, ...
    rangeVec)


nTarget = size(phiBounds,1);
objCount = 1;
for t = 1:nTarget
    facet_phi1 = phiBounds(t,1);
    facet_phi2 = phiBounds(t,2);
    objectHeight = heightVec(t);
    proposedRange = rangeVec(t);
    
    % handle foreground object
    deltaPhi = max([facet_phi1, facet_phi2])-min([facet_phi1, facet_phi2]); % angular extent of foreground object
    % halfWidth = tan(deltaPhi/2)*proposedRange;
    meanPhi = (facet_phi1 + facet_phi2)/2;
    
    edgeRange = proposedRange/cos(deltaPhi/2); % distance to edges of facet
    
    plane_point = [edgeRange*cos(facet_phi1), edgeRange*sin(facet_phi1), 0];
    plane_point2 = [edgeRange*cos(facet_phi2), edgeRange*sin(facet_phi2), 0];
    plane_span1 = plane_point2-plane_point; % points from plane_point to plain_point2;
    
    % plane_point = [-halfWidth, proposedRange, 0];
    % plane_span1 = [2*halfWidth,0,0];
    plane_span2 = [0,0,objectHeight];
    plane_normal = -[cos(meanPhi), sin(meanPhi), 0];
    
    proposedObjects{objCount,1} = 2;
    proposedObjects{objCount,2} = plane_point;
    proposedObjects{objCount,3} = plane_span1;
    proposedObjects{objCount,4} = plane_span2;
    proposedObjects{objCount,5} = plane_normal;
    objCount = objCount + 1;
end

end

