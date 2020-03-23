#!/bin/bash
MORTECH
for i in {1..80}
do
    cd defGeom/$i
    pimpleFoam
    cd ../..
done

