#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import pybind11
import subprocess

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

# A CMakeExtension needs a sourcedir instead of a file list.
# The name must be the _single_ output extension from the CMake build.
# If you need multiple extensions, see scikit-build.
class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))

        # required for auto-detection of auxiliary "native" libs
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        cfg = "Debug" if self.debug else "Release"
        # cfg = "Debug"

        # Set Python_EXECUTABLE instead if you use PYBIND11_FINDPYTHON
        # EXAMPLE_VERSION_INFO shows you how to pass a value into the C++ code
        # from Python.
        cmake_args = [
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={}".format(extdir),
            "-DPYTHON_EXECUTABLE={}".format(sys.executable),
            "-DCMAKE_BUILD_TYPE={}".format(cfg),  # not used on MSVC, but no harm
            "-DCMAKE_CXX_STANDARD=14"
        ]
        build_args = []

        # Enable required feature set in DDXspecial case for DDX:
        cmake_args += ["-DPYTHON=ON", "-DTESTS=OFF", "-DEXAMPLES=OFF",
                       "-DDDX_LIBRARY=OFF"]

        # Add Pybind11 info
        cmake_args += [f"-DPYBIND11_DIR={pybind11.get_cmake_dir()}"]

        # Set CMAKE_BUILD_PARALLEL_LEVEL to control the parallel build level
        # across all generators.
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ:
            # self.parallel is a Python 3 only way to set parallel jobs by hand
            # using -j in the build_ext call, not supported by pip or PyPA-build.
            if hasattr(self, "parallel") and self.parallel:
                # CMake 3.12+ only.
                build_args += ["-j{}".format(self.parallel)]

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        subprocess.check_call(
            ["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp
        )
        subprocess.check_call(
            ["cmake", "--build", "."] + build_args, cwd=self.build_temp
        )


def read_readme():
    with open("README.md") as fp:
        return "".join([line for line in fp if not line.startswith("<img")])


if not os.path.isfile("src/pyddx.cpp"):
    raise RuntimeError("Running setup.py is only supported "
                       "from top level of repository as './setup.py <command>'")

setup(
    ext_modules=[CMakeExtension("pyddx")],
    zip_safe=False,
    platforms=["Linux", "Mac OS-X"],
    cmdclass={"build_ext": CMakeBuild,}
)
