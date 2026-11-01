"""oligomet_deconv -- parallel charge-envelope deconvolution for oligonucleotide MS1 data.

Companion to the OligoMetProfiler R package: reads mzML/mzXML files streamed
scan-by-scan (bounded memory regardless of file size), detects chromatographic
features (ROIs), groups them into charge-state envelopes by neutral-mass
agreement across observed charge states, and optionally captures MS2 spectra
for a caller-supplied precursor watch-list. See cli.py for the entry point.
"""

__version__ = "0.1.0"
