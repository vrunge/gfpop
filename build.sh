#!/bin/bash
cd ..
set -o errexit
rm -rf gfpop-release
cp -r gfpop gfpop-release
PKG_TGZ=$(R CMD build gfpop-release|tee gfpop-release-build.txt|grep building|sed "s/.*\(gfpop.*.tar.gz\).*/\1/")
cat gfpop-release-build.txt
R --vanilla CMD INSTALL $PKG_TGZ
R --vanilla CMD check --as-cran $PKG_TGZ
