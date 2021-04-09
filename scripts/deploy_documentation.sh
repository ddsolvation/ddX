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
head=$(git rev-parse HEAD)
git config user.name "GitHub Actions Bot"
git config user.email "<>"

git fetch
git checkout -B gh-pages refs/remotes/origin/gh-pages

rm -rf dev
cp -a build/docs/html dev
git add dev
haschanges=$(git diff --cached dev)

if [ "$branch" != "main" ]; then
    echo "Skipping deployment as not on main."
    exit 0
fi
if [ -z "$haschanges" ]; then
    echo "Documentation did not change ... skipping deployment"
    exit 0
else
    git commit -m "Documentation build from $head"
    git push origin gh-pages
    git checkout $head
fi
