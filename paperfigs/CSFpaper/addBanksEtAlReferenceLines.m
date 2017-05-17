function addBanksEtAlReferenceLines(rParamsAllConds, plotLuminanceLine)
    sfIndex = 1;
    lumIndex = 1;
    rParams = rParamsAllConds{sfIndex,lumIndex};
    luminanceColors = [0 1 0; 0 0 1; 1 0 0; 0 0 0];
    luminanceColors = [0 0 0; 0 0 0; 0 0 0; 0 0 0];
    % Add unshifted version for reference
    banksFactor = 1;
    [A,B,C,D,E] = LoadDigitizedBanksFigure2;
    if (~rParams.oiParams.blur && ~rParams.mosaicParams.apertureBlur)
        plot(A(:,1),A(:,2),'k:','LineWidth',0.5);
        plot(A(:,1),A(:,2)*banksFactor,'r-','LineWidth',2);
    elseif (~rParams.oiParams.blur)
        plot(B(:,1),B(:,2),'k:','LineWidth',0.5);
        plot(B(:,1),B(:,2)*banksFactor,'r','LineWidth',2);
    else
        
        if (plotLuminanceLine(1))
            % 3.4 cd/m2
            plot(E(:,1),E(:,2)*banksFactor,'-','LineWidth',1.0, 'Color', squeeze(luminanceColors(1,:)));
        end
        
        if (plotLuminanceLine(2))
            % 34 cd/m2
            plot(D(:,1),D(:,2)*banksFactor,'-','LineWidth',1.0, 'Color', squeeze(luminanceColors(2,:)));
        end
        
        if (plotLuminanceLine(3))
            % 340 cd/m2
            plot(C(:,1),C(:,2)*banksFactor,'-','LineWidth',1.0, 'Color', squeeze(luminanceColors(3,:)));
        end
    end
end