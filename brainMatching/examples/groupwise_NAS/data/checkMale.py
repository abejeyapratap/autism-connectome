# Sanity check script to make sure all male subjects were written correctly
import os

dMales = []
sMales = []
# all male subject check
with open("desikan_male.txt", "r") as f:
    dMales = f.read().splitlines()

with open("schaefer_male.txt", "r") as f:
    sMales = f.read().splitlines()

print(len(dMales), len(sMales))

# ASD/TDC Check
with open("asd_desikan_male.txt", "r") as f:
    asdMales = f.read().splitlines()
with open("tdc_desikan_male.txt", "r") as f:
    tdcMales = f.read().splitlines()

print(len(tdcMales), len(asdMales))