function [Z, ncols, nrows, xllcorner, yllcorner, cellsize, nodata] = asciigridread( filename )
    %[Z, ncols, nrows, xllcorner, yllcorner, cellsize, nodata] = ...
    %asciigridread(FILENAME)
    %reads a grid from a file in Arc ASCII Grid format.
    %Z is a 2D array containing the data values.
    %The rest are all the items from the header.
    %NaN is assigned to elements of Z corresponding to nodata values
    %in the grid file.
    %Created 06/03/2007 by Bård Romstad (bard.romstad@geo.uio.no)
    %Last update 10/10/2007 


    % Open file.
    [fid, message] = fopen(filename,'r');
    if fid == -1
        error('%s: %s', filename, message);
    end

    % Read the 6-line header.
    [ncols, nrows, xllcorner, yllcorner, cellsize, nodata, msg] = readHeader(fid);
    if ~isempty(msg)
        fclose(fid);
        error(msg);
    end
% disp(nrows)
% disp(ncols)
    % Read the matrix of data values, putting the k-th row in the data
    % file into the k-th column of matrix Z.  Close file -- nothing left to
    % read after this.
    [Z, count] = fscanf(fid,'%f',[ncols,nrows]);
    fclose(fid);
% disp(count)
    % Replace each no-data value with NaN.
    Z(Z == nodata) = NaN;

    % Orient the data so that rows are parallel to the x-axis and columns
    % are parallel to the y-axis (for compatibility with MATLAT functions
    % like SURF and MESH).
    Z = Z';
end


%----------------------------------------------------------------------
function [ncols, nrows, xllcorner, yllcorner,...
                        cellsize, nodata, msg] = readHeader(fid)
    %function that reads the header of an asciigrid file and returns
    %the information in the header as variables
    itemNames = cell(1,6);
    value = cell(1,6);
    for k = 1:6
        line = fgetl(fid);
        [token,rem] = strtok(line);
        itemNames{k} = token;
        value{k} = str2double(rem);
    end

    ncols     = value{1};
    nrows     = value{2};
    xllcorner = value{3};
    yllcorner = value{4};
    cellsize  = value{5};
    nodata    = value{6};
    
    %check each item in the header
    msg = checkItemNames(itemNames);
end

%----------------------------------------------------------------------
function msg = checkItemNames(itemNames)
    %Function that checks each item in a list towards a predefined list
    msg = [];
    standardItems = {...
        'ncols',...
        'nrows',...
        'xllcorner',...
        'yllcorner',...
        'cellsize',...
        'nodata_value'};

    for k = 1:6
        if ~strcmpi(itemNames{k},standardItems{k})
            msg = sprintf('Unexpected item name ''%s'' in file header (line %d).',...
                          itemNames{k}, k);
            return;
        end
    end
end