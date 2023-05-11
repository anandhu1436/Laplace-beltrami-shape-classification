# Bilateral mesh denoising
The code uses libgl library

use command s to build and run the code

cmake -B build
cd build
make

To run a with a mesh from input folder

./bilateral ../input/Noisy.obj

use keybord to add and remove noise

N-add noise along normal
R-add randoem noise around vertex
B-one iteration of bilateral filter
S-save obj in denoise folder
W-save obj in noise folder
