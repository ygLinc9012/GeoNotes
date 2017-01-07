function quatMat = quat(vec)
    a = vec(1); 
    b = vec(2); 
    c = vec(3); 
    d = vec(4);
    quatMat =[ a -b -c -d; b  a -d  c; c  d  a -b; d -c  b  a];
end