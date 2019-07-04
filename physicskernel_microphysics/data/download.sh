#! /bin/bash -x

# please get from http://r-ccs-climate.riken.jp/members/yashiro/download/snapshot.microphysics.pe000003
#                 http://r-ccs-climate.riken.jp/members/yashiro/download/check.microphysics.pe000003

file=snapshot.microphysics.pe000003

rm -f ./$file
wget http://r-ccs-climate.riken.jp/members/yashiro/download/$file

echo "Checking md5sum:"
md5sum -c ${file}.md5 || exit 1

file=check.microphysics.pe000003

rm -f ./$file
wget http://r-ccs-climate.riken.jp/members/yashiro/download/$file

echo "Checking md5sum:"
md5sum -c ${file}.md5 || exit 1
