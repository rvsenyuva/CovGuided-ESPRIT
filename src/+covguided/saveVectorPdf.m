function saveVectorPdf(figOrAxes, filename)
%SAVE_VECTOR_PDF Export a MATLAB figure/axes as a **vector** PDF with embedded fonts.
%   save_vector_pdf(gcf, 'figure.pdf');
%   save_vector_pdf(gca, 'figure.pdf');
%
% - Forces 'painters' renderer (vector back-end).
% - Uses exportgraphics(..., 'ContentType','vector') to avoid rasterization.
% - Sets transparent background, suitable for IEEE journal submissions.
%
    if nargin < 1 || isempty(figOrAxes), figOrAxes = gcf; end
    if nargin < 2 || isempty(filename)
        error('save_vector_pdf:MissingFilename','Filename required');
    end

    % Resolve figure handle
    if ishghandle(figOrAxes, 'axes')
        h = ancestor(figOrAxes, 'figure');
    elseif ishghandle(figOrAxes, 'figure')
        h = figOrAxes;
    else
        error('save_vector_pdf:BadHandle','First argument must be a figure or axes handle.');
    end

    % Enforce vector-friendly settings
    set(h, 'Renderer', 'painters', 'Color','w', 'PaperPositionMode','auto', 'InvertHardcopy','off');

    try
        exportgraphics(h, filename, 'ContentType','vector', 'BackgroundColor','none');
        fprintf('Saved vector PDF to %s\n', filename);
    catch ME
        % Use a message ID and a format specifier (MATLAB requirement).
        % Also include the original identifier in the message body.
        warning('save_vector_pdf:ExportFailed', '[%s] %s', ME.identifier, ME.message);
    end
end
