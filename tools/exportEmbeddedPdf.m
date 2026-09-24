function exportEmbeddedPdf(fig, pdfPath, targetSize)
%EXPORTEMBEDDEDPDF Vector PDF with embedded fonts at an exact physical size.
%
% print -dpdf never embeds the base-14 fonts. exportgraphics does embed, but
% on this Retina Mac it writes all geometry under a 0.5 device transform
% (cm 0.125 * 4) while leaving font point sizes, line widths and marker
% sizes at their nominal values. Rendering the figure at 2x in every one of
% those quantities therefore produces the intended physical size with every
% glyph, stroke and marker at its intended size, and fonts embedded.
s = 2;
fig.Units = 'inches';
if nargin < 3 || isempty(targetSize)
    targetSize = fig.Position(3:4);
end
fig.Position(3:4) = targetSize * s;
objects = findall(fig);
% Read every value first, then write, so that changing a parent's font size
% cannot cascade into a child before the child is scaled itself.
props = {'FontSize', 'LineWidth', 'MarkerSize', 'ItemTokenSize'};
values = cell(numel(objects), numel(props));
for i = 1:numel(objects)
    for j = 1:numel(props)
        if isprop(objects(i), props{j})
            try values{i, j} = objects(i).(props{j}); catch, end
        end
    end
end
for i = 1:numel(objects)
    for j = 1:numel(props)
        if ~isempty(values{i, j}) && isnumeric(values{i, j})
            try objects(i).(props{j}) = values{i, j} * s; catch, end
        end
    end
end
% White full-canvas rectangle so exportgraphics' tight crop keeps the page.
annotation(fig, 'rectangle', [0, 0, 1, 1], 'Color', 'w', 'LineWidth', 0.01);
drawnow;
exportgraphics(fig, pdfPath, 'ContentType', 'vector');
% Guard against the device scale differing on another machine or release:
% the written page must match the requested size, not twice it.
fid = fopen(pdfPath, 'r'); bytes = fread(fid, '*uint8')'; fclose(fid);
tok = regexp(char(bytes), '/MediaBox\s*\[([^\]]+)\]', 'tokens', 'once');
b = sscanf(tok{1}, '%f'); got = [b(3)-b(1), b(4)-b(2)] / 72;
if any(abs(got - targetSize) > 0.05 * targetSize)
    warning('exportEmbeddedPdf:UnexpectedScale', ...
        ['Page is %.3f x %.3f in but %.3f x %.3f in was requested. The 0.5 ' ...
         'exportgraphics device scale this function compensates for does ' ...
         'not apply here; do not use the output.'], got, targetSize);
end
end
