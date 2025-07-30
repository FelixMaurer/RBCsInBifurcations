%--------------------------------------------------------------------------
% Script Name : inbounds.m
% Authors     : Felix Maurer
% Institution : Saarland University
% Email       : mail@felixmilanmaurer.com
% Date        : 2024
%
% Description :
%   This is a small helper function to check whether indices are within the
%   bounds of an array.
%
% Usage :
%
% Dependencies :
%
function flag = inbounds(A,x,y)
    sizeX = size(A,1);
    sizeY = size(A,2);
    if x > 0 && x <= sizeX && y > 0 && y <= sizeY
        flag = true;
    else
        flag = false;
    end
end