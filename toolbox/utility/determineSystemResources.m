function [numberOfCores, ramSizeGBytes, sizeofDoubleInBytes] = determineSystemResources()

   if (ismac)
        numberOfCores = feature('numCores');
        [s,q] = system('sysctl -n hw.memsize');
        ramSizeGBytes = str2double(q)/(1024*1024*1024);
   elseif (isunix)
        numberOfCores = feature('numCores');
        [s,q] = system(sprintf('free -h | gawk  ''/Mem:/{print $2}'''));
        if (s ~= 0)
            error('Cannot determine system resources because gawk is not installed. Please install using ''sudo apt install gawk''.');
        end
        ramSizeGBytes = str2double(strrep(q,'G',''));
   else
       error('No windows support');
   end

   aDouble = double(1.0);
   attr = whos('aDouble');
   sizeofDoubleInBytes = attr.bytes;
end
