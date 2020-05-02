function [jPhDrive] = f_photonDriveDef_threeByThreeClassifier(t,switchTimes,jPh0Vec)

jPhDrive = zeros(9,length(t));

for ii = 1:length(t)
    if t(ii) > switchTimes(1) && t(ii) <= switchTimes(2)
        
        jPhDrive(:,ii) = jPh0Vec(:);
        
    end
end