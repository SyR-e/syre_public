function downloadWithProgress(url, outFile)
    % DOWNLOADWITHPROGRESS Download a file with progress reporting.
    %
    %   downloadWithProgress(URL, OUTFILE)
    %   Streams URL into OUTFILE and prints progress in the Command Window.

    import matlab.net.http.*
    import matlab.net.http.io.*

    % HTTP GET request
    req = RequestMessage();

    % Stream directly to file (no large data in memory)
    consumer = FileConsumer(outFile);

    % Attach our custom progress monitor
    opts = HTTPOptions( ...
        'ProgressMonitorFcn', @MyProgressMonitor, ...
        'UseProgressMonitor', true);

    % Send request
    fprintf('Starting download from:\n  %s\n', url);
    resp = req.send(url, opts, consumer);

    % Basic error check
    if resp.StatusCode ~= StatusCode.OK
        error('Download failed: %s', char(resp.StatusLine));
    end

    fprintf('File saved to: %s\n', outFile);
end
