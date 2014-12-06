#!/bin/bash

cd bracelet
cd n3
~/repos/fftsvd/meshmaker N3.pdb ../radii.siz n3_stern_4_1.srf 1.4 2.0 4 1 1 1 .
cd ../n4
~/repos/fftsvd/meshmaker N4.pdb ../radii.siz n4_stern_4_1.srf 1.4 2.0 4 1 1 1 .
cd ../n5
~/repos/fftsvd/meshmaker N5.pdb ../radii.siz n5_stern_4_1.srf 1.4 2.0 4 1 1 1 .
cd ../n6
~/repos/fftsvd/meshmaker N6.pdb ../radii.siz n6_stern_4_1.srf 1.4 2.0 4 1 1 1 .
cd ../n7
~/repos/fftsvd/meshmaker N7.pdb ../radii.siz n7_stern_4_1.srf 1.4 2.0 4 1 1 1 .
cd ../n8
~/repos/fftsvd/meshmaker N8.pdb ../radii.siz n8_stern_4_1.srf 1.4 2.0 4 1 1 1 .

cd ../../rod/l5_q12
~/repos/fftsvd/meshmaker L5.pdb ../radii.siz L5_stern_4_1.srf 1.4 2.0 4 1 1 1 .
cd ../l5_q13
~/repos/fftsvd/meshmaker L5.pdb ../radii.siz L5_stern_4_1.srf 1.4 2.0 4 1 1 1 .
cd ../l5_q14
~/repos/fftsvd/meshmaker L5.pdb ../radii.siz L5_stern_4_1.srf 1.4 2.0 4 1 1 1 .
cd ../l5_q23
~/repos/fftsvd/meshmaker L5.pdb ../radii.siz L5_stern_4_1.srf 1.4 2.0 4 1 1 1 .
cd ../l6_q12
~/repos/fftsvd/meshmaker L6.pdb ../radii.siz L6_stern_4_1.srf 1.4 2.0 4 1 1 1 .
cd ../l6_q13
~/repos/fftsvd/meshmaker L6.pdb ../radii.siz L6_stern_4_1.srf 1.4 2.0 4 1 1 1 .
cd ../l6_q14
~/repos/fftsvd/meshmaker L6.pdb ../radii.siz L6_stern_4_1.srf 1.4 2.0 4 1 1 1 .
cd ../l6_q15
~/repos/fftsvd/meshmaker L6.pdb ../radii.siz L6_stern_4_1.srf 1.4 2.0 4 1 1 1 .


