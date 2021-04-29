function [xcoords, ycoords] = map_utaharray(fpath,electrode)
% create a channel Map file for simulated data (eMouse)

% here I know a priori what order my channels are in.  So I just manually
% make a list of channel indices (and give
% an index to dead channels too). chanMap(1) is the row in the raw binary
% file for the first channel. chanMap(1:2) = [33 34] in my case, which happen to
% be dead channels.

if nargin<2, electrode = 'utah96'; end % default

switch electrode
    case 'utah96'
        chanMap = 1:96;
        
        % the first thing Kilosort does is reorder the data with data = data(chanMap, :).
        % Now we declare which channels are "connected" in this normal ordering,
        % meaning not dead or used for non-ephys data
        
        connected = true(96, 1);
        
        % now we define the horizontal (x) and vertical (y) coordinates of these
        % 34 channels. For dead or nonephys channels the values won't matter. Again
        % I will take this information from the specifications of the probe. These
        % are in um here, but the absolute scaling doesn't really matter in the
        % algorithm.
        
        coords = [[3600, 3200];
            [3600, 2800];
            [3600, 2400];
            [3600, 2000];
            [3600, 1600];
            [3600, 1200];
            [3600, 800];
            [3600, 400];
            [3200, 3600];
            [3200, 3200];
            [3200, 2800];
            [3200, 2400];
            [3200, 2000];
            [3200, 1600];
            [3200, 1200];
            [3200, 800];
            [3200, 400];
            [3200, 0];
            [2800, 3600];
            [2800, 3200];
            [2800, 2800];
            [2800, 2400];
            [2800, 2000];
            [2800, 1600];
            [2800, 1200];
            [2800, 800];
            [2800, 400];
            [2800, 0];
            [2400, 3600];
            [2400, 3200];
            [2400, 2800];
            [2400, 2400];
            [2400, 2000];
            [2400, 1600];
            [2400, 1200];
            [2400, 800];
            [2400, 400];
            [2400, 0];
            [2000, 3600];
            [2000, 3200];
            [2000, 2800];
            [2000, 2400];
            [2000, 2000];
            [2000, 1600];
            [2000, 1200];
            [2000, 800];
            [2000, 400];
            [2000, 0];
            [1600, 3600];
            [1600, 3200];
            [1600, 2800];
            [1600, 2400];
            [1600, 2000];
            [1600, 1600];
            [1600, 1200];
            [1600, 800];
            [1600, 400];
            [1600, 0];
            [1200, 3600];
            [1200, 3200];
            [1200, 2800];
            [1200, 2400];
            [1200, 2000];
            [1200, 1600];
            [1200, 1200];
            [1200, 800];
            [1200, 400];
            [1200, 0];
            [800, 3600];
            [800, 3200];
            [800, 2800];
            [800, 2400];
            [800, 2000];
            [800, 1600];
            [800, 1200];
            [800, 800];
            [800, 400];
            [800, 0];
            [400, 3600];
            [400, 3200];
            [400, 2800];
            [400, 2400];
            [400, 2000];
            [400, 1600];
            [400, 1200];
            [400, 800];
            [400, 400];
            [400, 0];
            [0, 3200];
            [0, 2800];
            [0, 2400];
            [0, 2000];
            [0, 1600];
            [0, 1200];
            [0, 800];
            [0, 400]];
        
        xcoords = coords(:,1); xcoords = flipud(xcoords);
        ycoords = coords(:,2); ycoords = flipud(ycoords);
        
        % Often, multi-shank probes or tetrodes will be organized into groups of
        % channels that cannot possibly share spikes with the rest of the probe. This helps
        % the algorithm discard noisy templates shared across groups. In
        % this case, we set kcoords to indicate which group the channel belongs to.
        % In our case all channels are on the same shank in a single group so we
        % assign them all to group 1.
        
        kcoords = 1:96;
        
        % at this point in Kilosort we do data = data(connected, :), ycoords =
        % ycoords(connected), xcoords = xcoords(connected) and kcoords =
        % kcoords(connected) and no more channel map information is needed (in particular
        % no "adjacency graphs" like in KlustaKwik).
        % Now we can save our channel map for the eMouse.
        
        % would be good to also save the sampling frequency here
        fs = 30000;
        
        if nargin==1
            save(fullfile(fpath, 'chanMap.mat'), 'chanMap', 'connected', 'xcoords', 'ycoords', 'kcoords', 'fs');
        else
            xcoords = (xcoords/400)  + 1;
            ycoords = (ycoords/400)  + 1;
        end
    case 'utah2x48'
        chanMap = 1:96;
        
        % the first thing Kilosort does is reorder the data with data = data(chanMap, :).
        % Now we declare which channels are "connected" in this normal ordering,
        % meaning not dead or used for non-ephys data
        
        connected = true(96, 1);
        
        % now we define the horizontal (x) and vertical (y) coordinates of these
        % 34 channels. For dead or nonephys channels the values won't matter. Again
        % I will take this information from the specifications of the probe. These
        % are in um here, but the absolute scaling doesn't really matter in the
        % algorithm.
        
        coords = [[0, 0];
            [0, 400];
            [0, 800];
            [0, 1200];
            [0, 1600];
            [0, 2000];
            [400, 0];
            [400, 400];
            [400, 800];
            [400, 1200];
            [400, 1600];
            [400, 2000];
            [800, 0];
            [800, 400];
            [800, 800];
            [800, 1200];
            [800, 1600];
            [800, 2000];
            [1200, 0];
            [1200, 400];
            [1200, 800];
            [1200, 1200];
            [1200, 1600];
            [1200, 2000];
            [1600, 0];
            [1600, 400];
            [1600, 800];
            [1600, 1200];
            [1600, 1600];
            [1600, 2000];
            [2000, 0];
            [2000, 400];
            [2000, 800];
            [2000, 1200];
            [2000, 1600];
            [2000, 2000];
            [2400, 0];
            [2400, 400];
            [2400, 800];
            [2400, 1200];
            [2400, 1600];
            [2400, 2000];
            [2800, 0];
            [2800, 400];
            [2800, 800];
            [2800, 1200];
            [2800, 1600];
            [2800, 2000];
            [4000, 0];
            [4000, 400];
            [4000, 800];
            [4000, 1200];
            [4000, 1600];
            [4000, 2000];
            [4400, 0];
            [4400, 400];
            [4400, 800];
            [4400, 1200];
            [4400, 1600];
            [4400, 2000];
            [4800, 0];
            [4800, 400];
            [4800, 800];
            [4800, 1200];
            [4800, 1600];
            [4800, 2000];
            [5200, 0];
            [5200, 400];
            [5200, 800];
            [5200, 1200];
            [5200, 1600];
            [5200, 2000];
            [5600, 0];
            [5600, 400];
            [5600, 800];
            [5600, 1200];
            [5600, 1600];
            [5600, 2000];
            [6000, 0];
            [6000, 400];
            [6000, 800];
            [6000, 1200];
            [6000, 1600];
            [6000, 2000];
            [6400, 0];
            [6400, 400];
            [6400, 800];
            [6400, 1200];
            [6400, 1600];
            [6400, 2000];
            [6800, 0];
            [6800, 400];
            [6800, 800];
            [6800, 1200];
            [6800, 1600];
            [6800, 2000]];
        
        xcoords = coords(:,1);
        ycoords = coords(:,2);
        
        % Often, multi-shank probes or tetrodes will be organized into groups of
        % channels that cannot possibly share spikes with the rest of the probe. This helps
        % the algorithm discard noisy templates shared across groups. In
        % this case, we set kcoords to indicate which group the channel belongs to.
        % In our case all channels are on the same shank in a single group so we
        % assign them all to group 1.
        
        kcoords = 1:96;
        
        % at this point in Kilosort we do data = data(connected, :), ycoords =
        % ycoords(connected), xcoords = xcoords(connected) and kcoords =
        % kcoords(connected) and no more channel map information is needed (in particular
        % no "adjacency graphs" like in KlustaKwik).
        % Now we can save our channel map for the eMouse.
        
        % would be good to also save the sampling frequency here
        fs = 30000;
        
        if nargin==1
            save(fullfile(fpath, 'chanMap.mat'), 'chanMap', 'connected', 'xcoords', 'ycoords', 'kcoords', 'fs');
        else
            xcoords = (xcoords/400)  + 1;
            ycoords = (ycoords/400)  + 1;
        end
        
    case 'utah128'
        chanMap = 1:128;
        
        % the first thing Kilosort does is reorder the data with data = data(chanMap, :).
        % Now we declare which channels are "connected" in this normal ordering,
        % meaning not dead or used for non-ephys data
        
        connected = true(128, 1);
        
        % now we define the horizontal (x) and vertical (y) coordinates of these
        % 34 channels. For dead or nonephys channels the values won't matter. Again
        % I will take this information from the specifications of the probe. These
        % are in um here, but the absolute scaling doesn't really matter in the
        % algorithm.
        
        
        
        coords =   [[0    ,     400];
            [0     ,    800];
            [0     ,  1200];
            [0     ,   1600];
            [0     ,   2000];
            [0     ,   2400];
            [0     ,   2800];
            [400     ,   3200];
            [400     ,   3600];
            [400     ,   4000];
            [400     ,      0];
            [400     ,    400];
            [400     ,    800];
            [400     ,   1200];
            [400     ,   1600];
            [400     ,   2000];
            [800     ,   2400];
            [800     ,   2800];
            [800     ,   3200];
            [800     ,   3600];
            [800     ,   4000];
            [800     ,   4400];
            [800     ,      0];
            [800     ,    400];
            [800     ,    800];
            [800     ,   1200];
            [800     ,   1600];
            [1200     ,   2000];
            [1200     ,   2400];
            [1200     ,   2800];
            [1200     ,   3200];
            [1200     ,   3600];
            [1200     ,   4000];
            [1200     ,   4400];
            [1200     ,      0];
            [1200     ,    400];
            [1200     ,    800];
            [1200     ,   1200];
            [1200     ,  1600];
            [1600     ,   2000];
            [1600     ,   2400];
            [1600     ,   2800];
            [1600     ,  3200];
            [1600     ,   3600];
            [1600     ,   4000];
            [1600     ,   4400];
            [1600     ,      0];
            [1600     ,    400];
            [1600     ,    800];
            [1600     ,   1200];
            [1600     ,  1600];
            [2000     ,   2000];
            [2000     ,   2400];
            [2000     ,   2800];
            [2000     ,   3200];
            [2000     ,   3600];
            [2000     ,   4000];
            [2000     ,   4400];
            [2000     ,      0];
            [2000     ,   400];
            [2000     ,    800];
            [2000     ,  1200];
            [2000     ,  1600];
            [2400     ,   2000];
            [2400     ,   2400];
            [2400     ,   2800];
            [2400     ,   3200];
            [2400     ,   3600];
            [2400     ,   4000];
            [2400     ,   4400];
            [2400     ,      0];
            [2400     ,    400];
            [2400     ,    800];
            [2400     ,   1200];
            [2400     ,   1600];
            [2800     ,   2000];
            [2800     ,   2400];
            [2800     ,   2800];
            [2800     ,   3200];
            [2800     ,   3600];
            [2800     ,   4000];
            [2800     ,   4400];
            [2800     ,      0];
            [2800     ,    400];
            [2800     ,   800];
            [2800     ,   1200];
            [2800     ,   1600];
            [3200     ,   2000];
            [3200     ,   2400];
            [3200     ,   2800];
            [3200     ,   3200];
            [3200     ,   3600];
            [3200     ,   4000];
            [3200     ,   4400];
            [3200     ,      0];
            [3200     ,    400];
            [3200     ,    800];
            [3200     ,   1200];
            [3200     ,   1600];
            [3600     ,   2000];
            [3600     ,   2400];
            [3600     ,   2800];
            [3600     ,   3200];
            [3600     ,   3600];
            [3600     ,   4000];
            [3600     ,   4400];
            [3600     ,      0];
            [3600     ,   400];
            [3600     ,    800];
            [3600     ,   1200];
            [3600     ,   1600];
            [4000     ,   2000];
            [4000     ,   2400];
            [4000     ,   2800];
            [4000     ,   3200];
            [4000     ,   3600];
            [4000     ,   4000];
            [4000     ,   4400];
            [4000     ,    400];
            [4000     ,    800];
            [4400     ,   1200];
            [4400     ,  1600];
            [4400     ,  2000];
            [4400     ,  2400];
            [4400     ,  2800];
            [4400     ,  3200];
            [4400     ,  3600];
            [4400     ,  4000]];
        
        xcoords = coords(:,1); xcoords = flipud(xcoords);
        ycoords = coords(:,2); ycoords = flipud(ycoords);
        
        % Often, multi-shank probes or tetrodes will be organized into groups of
        % channels that cannot possibly share spikes with the rest of the probe. This helps
        % the algorithm discard noisy templates shared across groups. In
        % this case, we set kcoords to indicate which group the channel belongs to.
        % In our case all channels are on the same shank in a single group so we
        % assign them all to group 1.
        
        kcoords = 1:128;
        
        % at this point in Kilosort we do data = data(connected, :), ycoords =
        % ycoords(connected), xcoords = xcoords(connected) and kcoords =
        % kcoords(connected) and no more channel map information is needed (in particular
        % no "adjacency graphs" like in KlustaKwik).
        % Now we can save our channel map for the eMouse.
        
        % would be good to also save the sampling frequency here
        fs = 30000;
        
        if nargin==1
            save(fullfile(fpath, 'chanMap.mat'), 'chanMap', 'connected', 'xcoords', 'ycoords', 'kcoords', 'fs');
        else
            xcoords = (xcoords/400)  + 1;
            ycoords = (ycoords/400)  + 1;
        end
        
        %     case 'utah256' % for Ripple two 128ch Utah array recordings
        %
        %         chanMap = 1:256;
        %
        %         % the first thing Kilosort does is reorder the data with data = data(chanMap, :).
        %         % Now we declare which channels are "connected" in this normal ordering,
        %         % meaning not dead or used for non-ephys data
        %
        %         connected = true(256, 1);
        %
        %         % now we define the horizontal (x) and vertical (y) coordinates of these
        %         % 34 channels. For dead or nonephys channels the values won't matter. Again
        %         % I will take this information from the specifications of the probe. These
        %         % are in um here, but the absolute scaling doesn't really matter in the
        %         % algorithm.
        %
        %
        %
        %         coords =   [[0    ,     400];
        %             [0     ,    800];
        %             [0     ,  1200];
        %             [0     ,   1600];
        %             [0     ,   2000];
        %             [0     ,   2400];
        %             [0     ,   2800];
        %             [400     ,   3200];
        %             [400     ,   3600];
        %             [400     ,   4000];
        %             [400     ,      0];
        %             [400     ,    400];
        %             [400     ,    800];
        %             [400     ,   1200];
        %             [400     ,   1600];
        %             [400     ,   2000];
        %             [800     ,   2400];
        %             [800     ,   2800];
        %             [800     ,   3200];
        %             [800     ,   3600];
        %             [800     ,   4000];
        %             [800     ,   4400];
        %             [800     ,      0];
        %             [800     ,    400];
        %             [800     ,    800];
        %             [800     ,   1200];
        %             [800     ,   1600];
        %             [1200     ,   2000];
        %             [1200     ,   2400];
        %             [1200     ,   2800];
        %             [1200     ,   3200];
        %             [1200     ,   3600];
        %             [1200     ,   4000];
        %             [1200     ,   4400];
        %             [1200     ,      0];
        %             [1200     ,    400];
        %             [1200     ,    800];
        %             [1200     ,   1200];
        %             [1200     ,  1600];
        %             [1600     ,   2000];
        %             [1600     ,   2400];
        %             [1600     ,   2800];
        %             [1600     ,  3200];
        %             [1600     ,   3600];
        %             [1600     ,   4000];
        %             [1600     ,   4400];
        %             [1600     ,      0];
        %             [1600     ,    400];
        %             [1600     ,    800];
        %             [1600     ,   1200];
        %             [1600     ,  1600];
        %             [2000     ,   2000];
        %             [2000     ,   2400];
        %             [2000     ,   2800];
        %             [2000     ,   3200];
        %             [2000     ,   3600];
        %             [2000     ,   4000];
        %             [2000     ,   4400];
        %             [2000     ,      0];
        %             [2000     ,   400];
        %             [2000     ,    800];
        %             [2000     ,  1200];
        %             [2000     ,  1600];
        %             [2400     ,   2000];
        %             [2400     ,   2400];
        %             [2400     ,   2800];
        %             [2400     ,   3200];
        %             [2400     ,   3600];
        %             [2400     ,   4000];
        %             [2400     ,   4400];
        %             [2400     ,      0];
        %             [2400     ,    400];
        %             [2400     ,    800];
        %             [2400     ,   1200];
        %             [2400     ,   1600];
        %             [2800     ,   2000];
        %             [2800     ,   2400];
        %             [2800     ,   2800];
        %             [2800     ,   3200];
        %             [2800     ,   3600];
        %             [2800     ,   4000];
        %             [2800     ,   4400];
        %             [2800     ,      0];
        %             [2800     ,    400];
        %             [2800     ,   800];
        %             [2800     ,   1200];
        %             [2800     ,   1600];
        %             [3200     ,   2000];
        %             [3200     ,   2400];
        %             [3200     ,   2800];
        %             [3200     ,   3200];
        %             [3200     ,   3600];
        %             [3200     ,   4000];
        %             [3200     ,   4400];
        %             [3200     ,      0];
        %             [3200     ,    400];
        %             [3200     ,    800];
        %             [3200     ,   1200];
        %             [3200     ,   1600];
        %             [3600     ,   2000];
        %             [3600     ,   2400];
        %             [3600     ,   2800];
        %             [3600     ,   3200];
        %             [3600     ,   3600];
        %             [3600     ,   4000];
        %             [3600     ,   4400];
        %             [3600     ,      0];
        %             [3600     ,   400];
        %             [3600     ,    800];
        %             [3600     ,   1200];
        %             [3600     ,   1600];
        %             [4000     ,   2000];
        %             [4000     ,   2400];
        %             [4000     ,   2800];
        %             [4000     ,   3200];
        %             [4000     ,   3600];
        %             [4000     ,   4000];
        %             [4000     ,   4400];
        %             [4000     ,    400];
        %             [4000     ,    800];
        %             [4400     ,   1200];
        %             [4400     ,  1600];
        %             [4400     ,  2000];
        %             [4400     ,  2400];
        %             [4400     ,  2800];
        %             [4400     ,  3200];
        %             [4400     ,  3600];
        %             [4400     ,  4000];
        %             [0    ,     400];
        %             [0     ,    800];
        %             [0     ,  1200];
        %             [0     ,   1600];
        %             [0     ,   2000];
        %             [0     ,   2400];
        %             [0     ,   2800];
        %             [400     ,   3200];
        %             [400     ,   3600];
        %             [400     ,   4000];
        %             [400     ,      0];
        %             [400     ,    400];
        %             [400     ,    800];
        %             [400     ,   1200];
        %             [400     ,   1600];
        %             [400     ,   2000];
        %             [800     ,   2400];
        %             [800     ,   2800];
        %             [800     ,   3200];
        %             [800     ,   3600];
        %             [800     ,   4000];
        %             [800     ,   4400];
        %             [800     ,      0];
        %             [800     ,    400];
        %             [800     ,    800];
        %             [800     ,   1200];
        %             [800     ,   1600];
        %             [1200     ,   2000];
        %             [1200     ,   2400];
        %             [1200     ,   2800];
        %             [1200     ,   3200];
        %             [1200     ,   3600];
        %             [1200     ,   4000];
        %             [1200     ,   4400];
        %             [1200     ,      0];
        %             [1200     ,    400];
        %             [1200     ,    800];
        %             [1200     ,   1200];
        %             [1200     ,  1600];
        %             [1600     ,   2000];
        %             [1600     ,   2400];
        %             [1600     ,   2800];
        %             [1600     ,  3200];
        %             [1600     ,   3600];
        %             [1600     ,   4000];
        %             [1600     ,   4400];
        %             [1600     ,      0];
        %             [1600     ,    400];
        %             [1600     ,    800];
        %             [1600     ,   1200];
        %             [1600     ,  1600];
        %             [2000     ,   2000];
        %             [2000     ,   2400];
        %             [2000     ,   2800];
        %             [2000     ,   3200];
        %             [2000     ,   3600];
        %             [2000     ,   4000];
        %             [2000     ,   4400];
        %             [2000     ,      0];
        %             [2000     ,   400];
        %             [2000     ,    800];
        %             [2000     ,  1200];
        %             [2000     ,  1600];
        %             [2400     ,   2000];
        %             [2400     ,   2400];
        %             [2400     ,   2800];
        %             [2400     ,   3200];
        %             [2400     ,   3600];
        %             [2400     ,   4000];
        %             [2400     ,   4400];
        %             [2400     ,      0];
        %             [2400     ,    400];
        %             [2400     ,    800];
        %             [2400     ,   1200];
        %             [2400     ,   1600];
        %             [2800     ,   2000];
        %             [2800     ,   2400];
        %             [2800     ,   2800];
        %             [2800     ,   3200];
        %             [2800     ,   3600];
        %             [2800     ,   4000];
        %             [2800     ,   4400];
        %             [2800     ,      0];
        %             [2800     ,    400];
        %             [2800     ,   800];
        %             [2800     ,   1200];
        %             [2800     ,   1600];
        %             [3200     ,   2000];
        %             [3200     ,   2400];
        %             [3200     ,   2800];
        %             [3200     ,   3200];
        %             [3200     ,   3600];
        %             [3200     ,   4000];
        %             [3200     ,   4400];
        %             [3200     ,      0];
        %             [3200     ,    400];
        %             [3200     ,    800];
        %             [3200     ,   1200];
        %             [3200     ,   1600];
        %             [3600     ,   2000];
        %             [3600     ,   2400];
        %             [3600     ,   2800];
        %             [3600     ,   3200];
        %             [3600     ,   3600];
        %             [3600     ,   4000];
        %             [3600     ,   4400];
        %             [3600     ,      0];
        %             [3600     ,   400];
        %             [3600     ,    800];
        %             [3600     ,   1200];
        %             [3600     ,   1600];
        %             [4000     ,   2000];
        %             [4000     ,   2400];
        %             [4000     ,   2800];
        %             [4000     ,   3200];
        %             [4000     ,   3600];
        %             [4000     ,   4000];
        %             [4000     ,   4400];
        %             [4000     ,    400];
        %             [4000     ,    800];
        %             [4400     ,   1200];
        %             [4400     ,  1600];
        %             [4400     ,  2000];
        %             [4400     ,  2400];
        %             [4400     ,  2800];
        %             [4400     ,  3200];
        %             [4400     ,  3600];
        %             [4400     ,  4000]];
        %
        %
        %         xcoords = coords(:,1); xcoords = flipud(xcoords);
        %         ycoords = coords(:,2); ycoords = flipud(ycoords);
        %
        %         % Often, multi-shank probes or tetrodes will be organized into groups of
        %         % channels that cannot possibly share spikes with the rest of the probe. This helps
        %         % the algorithm discard noisy templates shared across groups. In
        %         % this case, we set kcoords to indicate which group the channel belongs to.
        %         % In our case all channels are on the same shank in a single group so we
        %         % assign them all to group 1.
        %
        %         kcoords = 1:256;
        %
        %         % at this point in Kilosort we do data = data(connected, :), ycoords =
        %         % ycoords(connected), xcoords = xcoords(connected) and kcoords =
        %         % kcoords(connected) and no more channel map information is needed (in particular
        %         % no "adjacency graphs" like in KlustaKwik).
        %         % Now we can save our channel map for the eMouse.
        %
        %         % would be good to also save the sampling frequency here
        %         fs = 30000;
        %
        %         if nargin==1
        %             save(fullfile(fpath, 'chanMap.mat'), 'chanMap', 'connected', 'xcoords', 'ycoords', 'kcoords', 'fs');
        %         else
        %             xcoords = (xcoords/400)  + 1;
        %             ycoords = (ycoords/400)  + 1;
        %         end
    case 'utah256' % for Ripple two 128ch Utah array recordings
        
        chanMap = 1:256;
        
        % the first thing Kilosort does is reorder the data with data = data(chanMap, :).
        % Now we declare which channels are "connected" in this normal ordering,
        % meaning not dead or used for non-ephys data
        
        connected = true(256, 1);
        
        % now we define the horizontal (x) and vertical (y) coordinates of these
        % 34 channels. For dead or nonephys channels the values won't matter. Again
        % I will take this information from the specifications of the probe. These
        % are in um here, but the absolute scaling doesn't really matter in the
        % algorithm.
        
        
        
        coords =   [[ 0, 0];
            [ 0 , 400];
            [ 0 , 800];
            [ 0 ,  1200];
            [ 0 ,  1600];
            [ 0 ,  2000];
            [ 0 ,  2400];
            [ 0 ,  2800];
            [ 0 ,  3200];
            [ 0 ,  3600];
            [ 0 ,  4000];
            [ 0 ,  4400];
            [ 0 ,  4800];
            [ 0 ,  5200];
            [ 0 ,  5600];
            [ 0 ,  6000];
            [ 0 ,  6400];
            [ 0 ,  6800];
            [ 0 ,  7200];
            [ 0 ,  7600];
            [ 0 ,  8000];
            [ 0 ,  8400];
            [ 0 ,  8800];
            [ 0 ,  9200];
            [ 0 ,  9600];
            [ 0 ,  10000];
            [ 0 ,  10400];
            [ 0 ,  10800];
            [ 0 ,  11200];
            [ 0 ,  11600];
            [ 0 ,  12000];
            [ 0 ,  12400];
            [ 0 ,  12800];
            [ 0 ,  13200];
            [ 0 ,  13600];
            [ 0 ,  14000];
            [ 0 ,  14400];
            [ 0 ,  14800];
            [ 0 ,  15200];
            [ 0 ,  15600];
            [ 0 ,  16000];
            [ 0 ,  16400];
            [ 0 ,  16800];
            [ 0 ,  17200];
            [ 0 ,  17600];
            [ 0 ,  18000];
            [ 0 ,  18400];
            [ 0 ,  18800];
            [ 0 ,  19200];
            [ 0 ,  19600];
            [ 0 ,  20000];
            [ 0 ,  20400];
            [ 0 ,  20800];
            [ 0 ,  21200];
            [ 0 ,  21600];
            [ 0 ,  22000];
            [ 0 ,  22400];
            [ 0 ,  22800];
            [ 0 ,  23200];
            [ 0 ,  23600];
            [ 0 ,  24000];
            [ 0 ,  24400];
            [ 0 ,  24800];
            [ 0 ,  25200];
            [ 0 ,  25600];
            [ 0 ,  26000];
            [ 0 ,  26400];
            [ 0 ,  26800];
            [ 0 ,  27200];
            [ 0 ,  27600];
            [ 0 ,  28000];
            [ 0 ,  28400];
            [ 0 ,  28800];
            [ 0 ,  29200];
            [ 0 ,  29600];
            [ 0 ,  30000];
            [ 0 ,  30400];
            [ 0 ,  30800];
            [ 0 ,  31200];
            [ 0 ,  31600];
            [ 0 ,  32000];
            [ 0 ,  32400];
            [ 0 ,  32800];
            [ 0 ,  33200];
            [ 0 ,  33600];
            [ 0 ,  34000];
            [ 0 ,  34400];
            [ 0 ,  34800];
            [ 0 ,  35200];
            [ 0 ,  35600];
            [ 0 ,  36000];
            [ 0 ,  36400];
            [ 0 ,  36800];
            [ 0 ,  37200];
            [ 0 ,  37600];
            [ 0 ,  38000];
            [ 0 ,  38400];
            [ 0 ,  38800];
            [ 0 ,  39200];
            [ 0 ,  39600];
            [ 0 ,  40000];
            [ 0 ,  40400];
            [ 0 ,  40800];
            [ 0 ,  41200];
            [ 0 ,  41600];
            [ 0 ,  42000];
            [ 0 ,  42400];
            [ 0 ,  42800];
            [ 0 ,  43200];
            [ 0 ,  43600];
            [ 0 ,  44000];
            [ 0 ,  44400];
            [ 0 ,  44800];
            [ 0 ,  45200];
            [ 0 ,  45600];
            [ 0 ,  46000];
            [ 0 ,  46400];
            [ 0 ,  46800];
            [ 0 ,  47200];
            [ 0 ,  47600];
            [ 0 ,  48000];
            [ 0 ,  48400];
            [ 0 ,  48800];
            [ 0 ,  49200];
            [ 0 ,  49600];
            [ 0 ,  50000];
            [ 0 ,  50400];
            [ 0 ,  50800];
            [ 0 ,  51200];
            [ 0 ,  51600];
            [ 0 ,  52000];
            [ 0 ,  52400];
            [ 0 ,  52800];
            [ 0 ,  53200];
            [ 0 ,  53600];
            [ 0 ,  54000];
            [ 0 ,  54400];
            [ 0 ,  54800];
            [ 0 ,  55200];
            [ 0 ,  55600];
            [ 0 ,  56000];
            [ 0 ,  56400];
            [ 0 ,  56800];
            [ 0 ,  57200];
            [ 0 ,  57600];
            [ 0 ,  58000];
            [ 0 ,  58400];
            [ 0 ,  58800];
            [ 0 ,  59200];
            [ 0 ,  59600];
            [ 0 ,  60000];
            [ 0 ,  60400];
            [ 0 ,  60800];
            [ 0 ,  61200];
            [ 0 ,  61600];
            [ 0 ,  62000];
            [ 0 ,  62400];
            [ 0 ,  62800];
            [ 0 ,  63200];
            [ 0 ,  63600];
            [ 0 ,  64000];
            [ 0 ,  64400];
            [ 0 ,  64800];
            [ 0 ,  65200];
            [ 0 ,  65600];
            [ 0 ,  66000];
            [ 0 ,  66400];
            [ 0 ,  66800];
            [ 0 ,  67200];
            [ 0 ,  67600];
            [ 0 ,  68000];
            [ 0 ,  68400];
            [ 0 ,  68800];
            [ 0 ,  69200];
            [ 0 ,  69600];
            [ 0 ,  70000];
            [ 0 ,  70400];
            [ 0 ,  70800];
            [ 0 ,  71200];
            [ 0 ,  71600];
            [ 0 ,  72000];
            [ 0 ,  72400];
            [ 0 ,  72800];
            [ 0 ,  73200];
            [ 0 ,  73600];
            [ 0 ,  74000];
            [ 0 ,  74400];
            [ 0 ,  74800];
            [ 0 ,  75200];
            [ 0 ,  75600];
            [ 0 ,  76000];
            [ 0 ,  76400];
            [ 0 ,  76800];
            [ 0 ,  77200];
            [ 0 ,  77600];
            [ 0 ,  78000];
            [ 0 ,  78400];
            [ 0 ,  78800];
            [ 0 ,  79200];
            [ 0 ,  79600];
            [ 0 ,  80000];
            [ 0 ,  80400];
            [ 0 ,  80800];
            [ 0 ,  81200];
            [ 0 ,  81600];
            [ 0 ,  82000];
            [ 0 ,  82400];
            [ 0 ,  82800];
            [ 0 ,  83200];
            [ 0 ,  83600];
            [ 0 ,  84000];
            [ 0 ,  84400];
            [ 0 ,  84800];
            [ 0 ,  85200];
            [ 0 ,  85600];
            [ 0 ,  86000];
            [ 0 ,  86400];
            [ 0 ,  86800];
            [ 0 ,  87200];
            [ 0 ,  87600];
            [ 0 ,  88000];
            [ 0 ,  88400];
            [ 0 ,  88800];
            [ 0 ,  89200];
            [ 0 ,  89600];
            [ 0 ,  90000];
            [ 0 ,  90400];
            [ 0 ,  90800];
            [ 0 ,  91200];
            [ 0 ,  91600];
            [ 0 ,  92000];
            [ 0 ,  92400];
            [ 0 ,  92800];
            [ 0 ,  93200];
            [ 0 ,  93600];
            [ 0 ,  94000];
            [ 0 ,  94400];
            [ 0 ,  94800];
            [ 0 ,  95200];
            [ 0 ,  95600];
            [ 0 ,  96000];
            [ 0 ,  96400];
            [ 0 ,  96800];
            [ 0 ,  97200];
            [ 0 ,  97600];
            [ 0 ,  98000];
            [ 0 ,  98400];
            [ 0 ,  98800];
            [ 0 ,  99200];
            [ 0 ,  99600];
            [ 0 ,  100000];
            [ 0 ,  100400];
            [ 0 ,  100800];
            [ 0 ,  101200];
            [ 0 ,  101600];
            [ 0 ,  102000]];
        
        
        xcoords = coords(:,1); xcoords = flipud(xcoords);
        ycoords = coords(:,2); ycoords = flipud(ycoords);
        
        % Often, multi-shank probes or tetrodes will be organized into groups of
        % channels that cannot possibly share spikes with the rest of the probe. This helps
        % the algorithm discard noisy templates shared across groups. In
        % this case, we set kcoords to indicate which group the channel belongs to.
        % In our case all channels are on the same shank in a single group so we
        % assign them all to group 1.
        
        kcoords = 1:256;
        
        % at this point in Kilosort we do data = data(connected, :), ycoords =
        % ycoords(connected), xcoords = xcoords(connected) and kcoords =
        % kcoords(connected) and no more channel map information is needed (in particular
        % no "adjacency graphs" like in KlustaKwik).
        % Now we can save our channel map for the eMouse.
        
        % would be good to also save the sampling frequency here
        fs = 30000;
        
        if nargin==1
            save(fullfile(fpath, 'chanMap.mat'), 'chanMap', 'connected', 'xcoords', 'ycoords', 'kcoords', 'fs');
        else
            xcoords = (xcoords/400)  + 1;
            ycoords = (ycoords/400)  + 1;
        end
end