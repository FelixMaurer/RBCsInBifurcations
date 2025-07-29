function flag = inbounds(A,x,y)
    sizeX = size(A,1);
    sizeY = size(A,2);
    if x > 0 && x <= sizeX && y > 0 && y <= sizeY
        flag = true;
    else
        flag = false;
    end
end