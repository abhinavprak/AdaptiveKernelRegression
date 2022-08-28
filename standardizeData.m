function[stdData, stdStat] = standardizeData(data, ref)
    ref = convertCharsToStrings(ref);
    if isstring(ref) 
        if ref == "self"
            location = mean(data);
            scale = std(data);
        end
    else
        location = ref.location;
        scale = ref.scale;
    end
    stdData = zeros(size(data));
    for i = 1:size(stdData,2)
        stdData(:,i) = (data(:,i) - location(i))/scale(i);
    end
    if nargout > 1
        stdStat.location = location;
        stdStat.scale = scale;
    end
end