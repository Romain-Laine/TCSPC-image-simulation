function  SaveAsOMETIFF(Data, Folder, dt)
% Adapted from code developped by Ian Munro "Munro, Ian" <i.munro@imperial.ac.uk>
% 2017-05-10

tic
disp('Saving as OME-tif...');
% verify that enough memory is allocated
bfCheckJavaMemory();
autoloadBioFormats = 1;

% load the Bio-Formats library into the MATLAB environment
status = bfCheckJavaPath(autoloadBioFormats);
assert(status, ['Missing Bio-Formats library. Either add loci_tools.jar '...
            'to the static Java path or add it to the Matlab path.']);

% initialize logging
loci.common.DebugTools.enableLogging('ERROR');
sizet = size(Data, 3);

Data_OME = reshape(Data, size(Data,1), size(Data,2), 1,1, size(Data,3));
% NB this line has been found to be crucial
java.lang.System.setProperty('javax.xml.transform.TransformerFactory', 'com.sun.org.apache.xalan.internal.xsltc.trax.TransformerFactoryImpl');

metadata = createMinimalOMEXMLMetadata(Data_OME);
modlo = loci.formats.CoreMetadata();

modlo.moduloT.type = loci.formats.FormatTools.LIFETIME;
modlo.moduloT.unit = 'ps';
% replace with 'Gated' if appropriate
modlo.moduloT.typeDescription = 'TCSPC';
modlo.moduloT.start = 0;

modlo.moduloT.step = dt;
modlo.moduloT.end = (sizet -1) * dt;

OMEXMLService = loci.formats.services.OMEXMLServiceImpl();
OMEXMLService.addModuloAlong(metadata,modlo,0);

% important to delete old versions before writing.
outputPath = [Folder, '.ome.tif'];
if exist(outputPath, 'file') == 2
    delete(outputPath);
end
bfsave(Data_OME, outputPath, 'metadata', metadata);
toc



end

