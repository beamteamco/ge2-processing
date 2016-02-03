#ge2-processing
A hexrd, numpy, scipy based collection of Python utilities to process GE2 files.

Currently spot IDing and cleaning is implemented.

To run:
`python smmoth_ge_data.py config.yml`

This will read the GE2 files specified in the config.yml file and detect centers of the spots.
Then it will write the spots to a new GE2 file named `synth_spots.ge2`.
