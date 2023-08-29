PyTrios
=======

PyTrios module to communicate with TriOS optical sensors over a serial port

G1 sensor range supported: 
-RAMSES SAM and SAMIP family of sensors (only spectrometer modules, inclination/pressure data not interpreted at present)
-MicroFlu sensors (e.g Chl, PC, CDOM)
-IPS boxes 

G2 sensor range supported/tested:
-RAMSES SAM with built in inclination/pressure

Additionally, a script providing a minimum command-line interface to trigger RAMSES readings and store to file, is provided. Example code on how to convert raw to calibrated readings is also included, this requires the instrument calibration files to be present.

This software is not an official TriOS product. For official TriOS software please visit http://www.trios.de/

Support from TriOS during implementation of the communication protocol is gratefully acknowledged.


No warranties, etc. You may choose to interpret the GPL license also as beerware. Enjoy!
