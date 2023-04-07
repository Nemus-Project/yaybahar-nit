# yaybahar-nit

Companion page for the paper *Efficient Simulation of the Yaybahar using a Modal Approach*.

This repository stores audio samples produced with the algorithms described into the paper. Each folder contains the audio files obtained with the three subsystems. The file names indicate the note played, while the suffix *_rauc* indicates that the bow force was increased to produce raucous sounds.

`BowedString` stores the audio samples obtained by running the bowed string algorithm. Thee are then fed into the spring as described into the paper.

`Spring` stores the audio samples obtained by running the spring algorithm, the file called *impulse* is a simple spring impulse response. The audio files are then fed into the membrane as described into the paper. Im

`Membrane` stores the audio samples obtained by running the membrane algorithm, the file called *impulse* was obtained by feeding the spring impulse response into the membrane. This folder stores the final output of the instrument.

