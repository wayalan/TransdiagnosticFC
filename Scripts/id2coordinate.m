function coordinates = id2coordinate(id, numROI)
% convert one-dimension id to 2 dimension coordinates
t = 0;
coordinatesAll = zeros((numROI*numROI-numROI)/2,2);
for i = 1 : numROI-1
    for j = i+1 : numROI
        t = t + 1;
        coordinatesAll(t,:) = [i,j];
    end
end
coordinates = coordinatesAll(id,:);