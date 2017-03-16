function dirname = paramsToInstanceDirName(obj,params)
% dirname = paramsToInstanceDirName(obj,params)
% 
% Generate a directory names that captures the basic LMPlane instance generation
% parameters.

if (~strcmp(params.type,'Instance'))
    error('Incorrect parameter type passed');
end

    switch(params.instanceType)
        case 'LMPlane'
            theLMPlaneInstanceName = sprintf('LM_angles%0.0f_to_%0.0f_N%0.0f_contrasts%0.5f_to_%0.3f_N%0.0f_trialsN%0.0f',...
                params.startAngle,...
                params.deltaAngle,...
                params.nAngles,...
                params.lowContrast, ...
                params.highContrast, ...
                params.nContrastsPerDirection, ...
                params.trialsNum ...
                );
            dirname = theLMPlaneInstanceName;
            
        case 'LMSPlane'
            theLMSPlaneInstanceName = sprintf('LMS_azimuth%0.0f_to_%0.0f_N%0.0f_elevation%0.0f_to_%0.0f_N%0.0fcontrasts%0.5f_to_%0.3f_N%0.0f_trialsN%0.0f',...
                params.startAzimuthAngle,...
                params.deltaAzimuthAngle,...
                params.nAzimuthAngles,...
                params.startElevationAngle,...
                params.deltaElevationAngle,...
                params.nElevationAngles,...
                params.lowContrast, ...
                params.highContrast, ...
                params.nContrastsPerDirection, ...
                params.trialsNum ...
                );
            dirname = theLMSPlaneInstanceName;
            
        case 'contrasts'
            theContrastInstanceName = sprintf('C_contrasts%0.5f_to_%0.3f_N%0.0f_trialsN%0.0f',...
                params.lowContrast, ...
                params.highContrast, ...
                params.nContrastsPerDirection, ...
                params.trialsNum ...
                );
            dirname = theContrastInstanceName;
        
        case 'pedestalIncrements'
            thePedestalInstanceName = sprintf('pedestal_contrasts%0.5f_to_%0.3f_N%0.0f_trialsN%0.0f',...
                params.lowContrast, ...
                params.highContrast, ...
                params.nContrastsPerDirection, ...
                params.trialsNum ...
                );
            dirname = thePedestalInstanceName;
            
        otherwise
            error('Unknown instance type passed: ''%s''.', params.instanceType);  
    end
end




