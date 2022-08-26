proposedObject = makeProposedObject_multiTarget(phiBoundsHat, hhat, ...
    rhat);


figure(4)
plotObjects(proposedObject,'blue')
if ~isempty(objects)
    plotObjects(objects(1:3,:),'red')
end
hold off

fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperPosition = 4*[-.3 -.1 3 3];
fig.Units = 'centimeters';
fig.PaperSize=4*[2.45 2.8];
print(fig,[dirpath 'fgSceneEstimate_' num2str(posNum) ], '-dpdf','-r200');

