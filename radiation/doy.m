function res = doy(date)
    DV  = datevec(date);  % [N x 6] array
    DV  = DV(:, 1:3);   % [N x 3] array, no time
    DV2 = DV;
    DV2(:, 2:3) = 0;    % [N x 3], day before 01.Jan
    res = datenum(DV) - datenum(DV2);
end