#!/bin/bash

cd bracelet
cd n3
~/repos/fftsvd/meshmaker N3.pdb ../radii.siz n3_6.srf 1.4 2.0 6 1 1 0 .
cd ../n4
~/repos/fftsvd/meshmaker N4.pdb ../radii.siz n4_6.srf 1.4 2.0 6 1 1 0 .
cd ../n5
~/repos/fftsvd/meshmaker N5.pdb ../radii.siz n5_6.srf 1.4 2.0 6 1 1 0 .
cd ../n6
~/repos/fftsvd/meshmaker N6.pdb ../radii.siz n6_6.srf 1.4 2.0 6 1 1 0 .
cd ../n7
~/repos/fftsvd/meshmaker N7.pdb ../radii.siz n7_6.srf 1.4 2.0 6 1 1 0 .
cd ../n8
~/repos/fftsvd/meshmaker N8.pdb ../radii.siz n8_6.srf 1.4 2.0 6 1 1 0 .

cd ../../rod/l5_q12
~/repos/fftsvd/meshmaker L5.pdb ../radii.siz L5_6.srf 1.4 2.0 6 1 1 0 .
cd ../l5_q13
~/repos/fftsvd/meshmaker L5.pdb ../radii.siz L5_6.srf 1.4 2.0 6 1 1 0 .
cd ../l5_q14
~/repos/fftsvd/meshmaker L5.pdb ../radii.siz L5_6.srf 1.4 2.0 6 1 1 0 .
cd ../l5_q23
~/repos/fftsvd/meshmaker L5.pdb ../radii.siz L5_6.srf 1.4 2.0 6 1 1 0 .
cd ../l6_q12
~/repos/fftsvd/meshmaker L6.pdb ../radii.siz L6_6.srf 1.4 2.0 6 1 1 0 .
cd ../l6_q13
~/repos/fftsvd/meshmaker L6.pdb ../radii.siz L6_6.srf 1.4 2.0 6 1 1 0 .
cd ../l6_q14
~/repos/fftsvd/meshmaker L6.pdb ../radii.siz L6_6.srf 1.4 2.0 6 1 1 0 .
cd ../l6_q15
~/repos/fftsvd/meshmaker L6.pdb ../radii.siz L6_6.srf 1.4 2.0 6 1 1 0 .


