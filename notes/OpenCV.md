## Schleiren notes  
OK, let's use 3 channels to store the edgefiltered pulse simulation and the etalon enhanced simulaitons
Base that output on the result of the various python work you've done lately
Use a temporary buffer, if that seems better than a fixed malloc.
instead of 3 single channels, use 1 3channel temporary buffer and fill it in the `pulsearray[f].fillrow_uint8c3` call
Then we need to come up with an optical scheme for doing the schleiren thing spectrally

This feels like it uses 3 monochrome cameras... sounds crazy.
For sake of PET though, schleiren might be the way to go for the imaging aspect.
