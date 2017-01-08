function A = triangulation2adjacency(face)
    f = double(face);

    A = sparse([f(:,1); f(:,1); f(:,2); f(:,2); f(:,3); f(:,3)], ...
               [f(:,2); f(:,3); f(:,1); f(:,3); f(:,1); f(:,2)], ...
               1.0);
    % avoid double links
    A = double(A>0);
end