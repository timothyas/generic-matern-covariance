# Source code

## Read in with `open_smoothdataset`

For reading in smooth output where all samples are written to a single file (now
deprecated) see [smooth_store.py](smooth_store.py).
This uses `xmitgcm.open_mdsdataset`, and does not perform well with LLC data.

For reading in samples where each sample is written with an independent 10-digit
"iteration number", see [new_smooth_store.py](new_smooth_store.py).
This uses `xmitgcm.llcreader` when `geometry='llc'`, but `xmitgcm.open_mdsdataset`
otherwise.

See [llcreader_with_smooth.ipynb](llcreader_with_smooth.ipynb) for some chunking
and performance notes.
Note that by default, this chunks `face` by 3 rather than 1, has
~4.2 tasks per chunk, `mds_store` does 5 ... but `llcreader` can instantly
(?) chunk all z levels at once, leading to far fewer chunks. Can't do this
with `mds_store`...
