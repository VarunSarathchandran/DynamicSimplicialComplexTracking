function exportFigurePdf(figureHandle, pdfPath)
%EXPORTFIGUREPDF Export a vector PDF at the figure's physical dimensions.

% Setting PaperSize explicitly avoids a macOS/Retina exportgraphics scaling
% issue that can halve both the requested page size and the rendered fonts.
validateattributes(figureHandle, {'matlab.ui.Figure'}, {'scalar'});
assert(ischar(pdfPath) || (isstring(pdfPath) && isscalar(pdfPath)), ...
    'dynsc:InvalidFigurePath', ...
    'pdfPath must be a character vector or scalar string.');

originalUnits = figureHandle.Units;
figureHandle.Units = 'inches';
paperSize = figureHandle.Position(3:4);
figureHandle.Units = originalUnits;

set(figureHandle, ...
    'PaperUnits', 'inches', ...
    'PaperPosition', [0, 0, paperSize], ...
    'PaperSize', paperSize);
print(figureHandle, char(pdfPath), '-dpdf', '-vector');
end
