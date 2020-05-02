function [jPhDrive] = f_photonDriveDef(t,switchTimes,jPh0)

jPhDrive = zeros(size(t));

for ii = 1:length(jPhDrive)
    for jj = 1:length(switchTimes)-1
        
        if t(ii) > switchTimes(jj) && t(ii) <= switchTimes(jj+1)
            jPhDrive(ii) = jPh0(jj);
        end
        
    end
end