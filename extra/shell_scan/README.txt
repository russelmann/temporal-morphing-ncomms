Description

Reading 3D scanning data for shells and measuring error in shape replication. This does not contain
raw scanning data. In order to run processing of raw data, download it from data repository and
place contents of the folder <shell scan> in the folder <data> here. Note that some of the
reconstructed marker locations are manually removed in cases when they are falsely detected (when,
in fact, there is no marker in the actual scanning data).

Folder <data>: Preprocessed 3D scanning data. Marker locations from the scans and from simulation.

File <scan_preprocess.m>: Preprocess raw 3D scanning data. Extract marker locations.

File <scan_register.m>: Register markers from scanning data to locations from simulation.

File <scan_register_weave.m>: Register markers from scanning data to locations from simulation.
Custom code for shell "weave" (Figure 4d).
