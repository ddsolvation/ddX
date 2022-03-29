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
rm -f .gitignore
cp -a build/docs/html dev
cd dev
git add .

if [ "$branch" != "main" ]; then
    echo "Skipping deployment as not on main."
    exit 0
fi
if git status | grep -q 'up to date'; then
    echo "Documentation did not update ... skipping deployment"
    exit 0
else
    git commit -m "Documentation build from $head"
    git push -f origin gh-pages
    git checkout $head
fi
