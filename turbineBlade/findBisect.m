function bisector = findBisect(p1, p2, p3)
    vec1 = p2 - p1;
    vec2 = p3 - p1;
    vec1 = vec1 / norm(vec1);
    vec2 = vec2/ norm(vec2);
    
    bisector = vec1 + vec2;
    bisector = bisector' / norm(bisector);
    
    cross_prod = cross([vec1], [vec2]);
    
    if cross_prod(1) > 0
        bisector = -bisector;
    end
    bisector = bisector';
end
