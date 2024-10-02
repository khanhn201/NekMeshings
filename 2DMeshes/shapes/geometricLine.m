function lines = geometricLine(N, R_start, R_end, multiplier, height)
    lines = [];
    idx = 0;
    series = geometricSeries(N + 1, R_start, R_end, multiplier);
    for i = 1:N
        idx += 1;
        lines(idx, 1, 1) = series(i);
        lines(idx, 1, 2) = height;
        lines(idx, 2, 1) = series(i+1);
        lines(idx, 2, 2) = height;
    end

end

function series = geometricSeries(N, R_start, R_end, multiplier)
    series = zeros(N, 1);
    series(1) = R_start;

    total_length = R_end - R_start;
    ratio = (multiplier^N - 1) / (multiplier - 1);
    spacing = total_length / ratio;
    for i = 2:N
        series(i) = series(i-1) + spacing;
        spacing = spacing*multiplier;
    end
end
