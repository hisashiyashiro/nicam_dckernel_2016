#! /bin/bash -x

# please get from http://r-ccs-climate.riken.jp/members/yashiro/download/snapshot.radiation.pe000003
#                 http://r-ccs-climate.riken.jp/members/yashiro/download/check.radiation.pe000003

file=snapshot.radiation.pe000003

rm -f ./$file
wget http://r-ccs-climate.riken.jp/members/yashiro/download/$file

echo "Checking md5sum:"
md5sum -c ${file}.md5 || exit 1

file=check.radiation.pe000003

rm -f ./$file
wget http://r-ccs-climate.riken.jp/members/yashiro/download/$file

echo "Checking md5sum:"
md5sum -c ${file}.md5 || exit 1
