#!/bin/bash
rm -rf output/lin.c-lin.bak
cp -r output/lin.c-lin output/lin.c-lin.bak
./quarkEigenLin.exe output/lin.c-lin/c-lin.cfg
