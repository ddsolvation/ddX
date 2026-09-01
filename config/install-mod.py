#!/usr/bin/env python3

from os import environ, listdir, makedirs
from os.path import exists, isdir, join
from shutil import copy
from sys import argv

build_dir = environ["MESON_BUILD_ROOT"]
if "MESON_INSTALL_DESTDIR_PREFIX" in environ:
    install_dir = environ["MESON_INSTALL_DESTDIR_PREFIX"]
else:
    install_dir = environ["MESON_INSTALL_PREFIX"]

include_dir = argv[1] if len(argv) > 1 else "include"
module_dir = join(install_dir, include_dir)

modules = []
for directory in listdir(build_dir):
    candidate = join(build_dir, directory)
    if isdir(candidate):
        for filename in listdir(candidate):
            if filename.endswith(".mod"):
                modules.append(join(candidate, filename))

if not exists(module_dir):
    makedirs(module_dir)

for module in modules:
    print("Installing", module, "to", module_dir)
    copy(module, module_dir)
