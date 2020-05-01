function [jPh0Vec] = f_createInputImageDrive(letter,defect,jPh0)

% see Nature 521 61 (2015)
% letter is 'z', 'v', or 'n'
% defect is integer 1-9, flips that pixel; defect = 0 gives ideal image
% for 3 x 3 image, indexing goes down the columns, like Fig. 2(a) of the
% Strukov paper

switch letter
    case 'z'
        jPh0Vec = [1 0 0 1 1 1 0 0 1];
    case 'v'
        jPh0Vec = [1 1 0 0 0 1 1 1 0];
    case 'n'
        jPh0Vec = [0 1 1 1 0 0 0 1 1];
end

if defect > 0
    if jPh0Vec(defect) == 0; jPh0Vec(defect) = 1; elseif jPh0Vec(defect) == 1; jPh0Vec(defect) = 0; end
end

jPh0Vec = jPh0*jPh0Vec;