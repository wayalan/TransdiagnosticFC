function coordinates = id2coordinate_dia(id, numROI)
% convert one-dimension id to 2 dimension coordinates
t = 0;
coordinatesAll = zeros(((numROI*numROI-numROI)/2)+numROI,2);
for i = 1 : numROI
    for j = i : numROI
        t = t + 1;
        coordinatesAll(t,:) = [i,j];
    end
end
coordinates = coordinatesAll(id,:);