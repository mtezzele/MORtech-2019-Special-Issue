#!/bin/bash
MORTECH
for i in {1..170}
do
    cd defGeom/$i
    pimpleFoam
    cd ../..
done

