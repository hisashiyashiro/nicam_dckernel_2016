#! /bin/bash -x

# please get from http://scale.aics.riken.jp/yashiro/download/snapshot.radiation.pe000003
#                 http://scale.aics.riken.jp/yashiro/download/check.radiation.pe000003

rm -f ./snapshot.radiation.pe000003
rm -f ./check.radiation.pe000003
wget http://scale.aics.riken.jp/yashiro/download/snapshot.radiation.pe000003
wget http://scale.aics.riken.jp/yashiro/download/check.radiation.pe000003
