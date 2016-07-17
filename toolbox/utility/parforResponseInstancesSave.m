function parforResponseInstancesSave(theStimData,filename)

save(filename, 'theStimData','-v7.3');
fprintf('Saved response instances to ''%s''.\n', filename);

end
