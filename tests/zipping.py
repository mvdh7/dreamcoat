import tarfile
import os

cpath = "tests/figures/"
files = [f for f in os.listdir(cpath) if f.endswith(".png")]

archive = "tests/data/archive.tar.gz"

with tarfile.open(archive, "w:xz") as tar:
    for file in files:
        tar.add(cpath + file)  # can also be dir (added recursively), symlink, etc


#%%
with tarfile.open("tests/figures/0.tar.gz", "w:xz") as tar:
    tar.add(cpath + files[0])

#%%

import bz2, os, sys

filename_in = "tests/figures/0_surphys_current_speed.png"
filename_out = "tests/figures/0_surphys_current_speed.bz2"

with open(filename_in, mode="rb") as fin, bz2.open(filename_out, "wb") as fout:
    fout.write(fin.read())

print(f"Uncompressed size: {os.stat(filename_in).st_size}")
# Uncompressed size: 1000000
print(f"Compressed size: {os.stat(filename_out).st_size}")
# Compressed size: 48
