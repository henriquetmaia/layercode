%useful for assessing accuracy, given files are formated as 
% decimal_base_code_viewpoint.png
function correctEncoding = decodeFilename( filename )
    [~, name, ~] = fileparts(filename);
    splits = strsplit(name, '_');
    decimalStr = splits(1);
    decimalNum = str2double(decimalStr);
    correctEncoding = de2bi(decimalNum,'left-msb');
end

