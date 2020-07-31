function c = myCal(rawdata, gain, offst)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


c = (rawdata .* gain) + offst;
end

