function visualizeResponses(coneMosaic, coneExcitations, figNo)
% Visualize the cone mosaic responses
%
% Syntax:
%   visualizeResponses(coneMosaic, coneExcitations, figNo)
%
% Description:
%    Visualize the cone mosaic responses on the provided figure, given the
%    cone mosaic and cone excitations.
%
% Inputs:
%    coneMosaic      - Object. A cone mosaic object.
%    coneExcitations - Matrix. A 3D matrix of cone excitations with a 2D
%                      matrix of cone excitations for each L, M, & S cones.
%    figNo           - Numeric. The figure number.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

    hFig = figure(figNo);
    clf;
    set(hFig, 'Position', [10 10 1500 950], 'Color', [1 1 1], ...
        'Name', 'Mosaic and 3 response instances');
    axHandle = subplot(2, 3, 2);
    coneMosaic.visualizeGrid('axesHandle', axHandle, ...
        'ticksInVisualDegs', true);
    xlabel('space (degs)', 'FontWeight', 'bold');
    ylabel('space (degs)', 'FontWeight', 'bold');
    title('cone mosaic');

    for k = 1:3
        axHandle = subplot(2, 3, 3 + k);
        if (k == 3)
            coneMosaic.renderActivationMap(...
                axHandle, squeeze(coneExcitations(k, :, :)), ...
                'mapType', 'modulated disks', ...
                'signalRange', [0 max(coneExcitations(:))], ...
                'showColorBar', true, ...
                'labelColorBarTicks', true, ...
                'titleForColorBar', sprintf('R*/cone/%2.0fms', ...
                coneMosaic.integrationTime * 1000));
        else
            coneMosaic.renderActivationMap(...
                axHandle, squeeze(coneExcitations(k, :, :)), ...
                'mapType', 'modulated disks', ...
                'signalRange', [0 max(coneExcitations(:))], ...
                'showColorBar', ~true, ...
                'labelColorBarTicks', ~true);
        end
        if (k > 1)
        	set(gca, 'YTickLabel', {});
            ylabel('');
        end
        xlabel('space (degs)', 'FontWeight', 'bold');
        title(sprintf('response instance #%d', k));
    end
end