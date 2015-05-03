clear,clc;

load('data\D.mat');
ComparisonMethods = D(1).methods;
nrCompMethods = size(ComparisonMethods,2);
maxarray = zeros(1,nrCompMethods);
minarray = ones(1,nrCompMethods)*inf;
maxOutliers = zeros(1,nrCompMethods);
minOutliers = zeros(1,nrCompMethods);
for j = 1:nrCompMethods
    for i = 1:size(D,2)
        if isinf(D(i).results(j)) 
            D(i).results(j) = 0;
        end
        if(D(i).results(j)>maxarray(j))
            maxarray(j) = D(i).results(j);
            maxOutliers(j) = i;
        end
        if(D(i).results(j)<minarray(j))
            minarray(j) = D(i).results(j);
            minOutliers(j) = i;
        end
    end
end

ComparisonMethods
minOutliers
maxOutliers