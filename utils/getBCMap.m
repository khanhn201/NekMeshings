function bcMap = getBCMap()
    ids   = [    1,     2,     3,     4,     5,     6,     7];
    names = {'W  ', 'v  ', 'int', 'o  ', 's  ', 'SYM', 'P  '};
    namesScalar = {'I  ', 't  ', 'int', 'O  ', 's  ', 'SYM', 'P  '};
    bcMap.id2str = containers.Map(ids, names);
    bcMap.str2id = containers.Map(names, ids);
    bcMap.id2strS = containers.Map(ids, namesScalar);
    bcMap.str2idS = containers.Map(namesScalar, ids);
end
