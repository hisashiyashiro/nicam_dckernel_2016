#! /bin/bash -x

# please get from http://r-ccs-climate.riken.jp/members/yashiro/download/snapshot.dc_vi_rhow_solver.pe000000

file=snapshot.dc_vi_rhow_solver.pe000000

rm -f ./$file
wget http://r-ccs-climate.riken.jp/members/yashiro/download/$file

echo "Checking md5sum:"
md5sum -c ${file}.md5 || exit 1
