The SEGY reader has been written by Aqeel Al-Naser <aqeel.al-naser@cs.manchester.ac.uk>.
Please see attached files. Please note the following assumptions and limitations for the SEGY reader:

1. No support for SEGY with 4-byte IBM floating-point sample format. It needs special treatment.
   However, it is not common. The other supported formats are:
     - 4-byte IEEE floating-point
     - 4-byte, two's complement integer
     - 2-byte, two's complement integer
     - 1-byte, two's complement integer
2. Current version assumes the followings in regards to the handled SEGY file:
     - no extended textual header records
     - all traces have same number of samples
     - all in-lines have same number of x-lines
3. In this version, dimension of SEGY file is calculated in an expensive way as it is not included in
   the header of the file; it has to go through all samples in the file. There might be an efficient way though.
4. The following default values are employed:
     - spacing: (1,1,1)
     - transformation: identity
     - time step: -1.0f
     - modality: MODALITY_UNKNOWN
     - slice order: +z
5. readBrick() function is not implemented
6. Exception handling needs integration with Voreen environment as is the case with RawVolumeReader.
7. As in RawVolumeReader, float samples are converted into Int16! This gives a better rendering results but not sure why!
8. Progress per upload step is not updated.
