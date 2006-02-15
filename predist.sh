#!/bin/sh -x

# build article
(
cd article
make
)

rm -f predist.sh
