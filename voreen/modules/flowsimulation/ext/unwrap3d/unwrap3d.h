//This program is written by Hussein Abdul-Rahman and Munther Gdeisat to program the three-dimensional phase unwrapper
//entitled "Fast three-dimensional phase-unwrapping algorithm based on sorting by
//reliability following a noncontinuous path"
//by  Hussein Abdul-Rahman, Munther A. Gdeisat, David R. Burton, and Michael J. Lalor,
//published in the Applied Optics, Vol. 46, No. 26, pp. 6623-6635, 2007.
//This program was written on 29th August 2007

//The wrapped phase volume is floating point data type. Also, the unwrapped phase volume is foloating point
//read the data from the file frame by frame
//The mask is byte data type.
//When the mask is 1 (true)  this means that the voxel is valid
//When the mask is 0 (false) this means that the voxel is invalid (noisy or corrupted voxel)
void unwrap3d(const float* WrappedVolume, float* UnwrappedVolume, const unsigned char* mask, int volume_width, int volume_height, int volume_depth);