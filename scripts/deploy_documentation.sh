#!/bin/bash -e

THISDIR=$(dirname "${BASH_SOURCE[0]}")
cd "$THISDIR/.."

if [ -d build ]; then
	echo "Delete build dir before running this script" >&2
	exit 1
fi

mkdir build
pushd build
cmake ..
make docs
popd

branch=$(git branch --show-current)
if [ "$branch" != "main" ]; then
	echo "Skipping deployment as not on main."
	exit 0
fi

head=$(git rev-parse HEAD)
git checkout gh-pages
rm -rf dev
cp -a build/docs/html dev
git add dev
git commit -m "Documentation build from $head" dev
git push -f
git checkout $head
