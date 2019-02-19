#!/bin/bash

find . -type f -printf 'touch --no-create -d "%t" "%p"\n' > ${IEMIC_BUILD_DIR}/restore_timestamps.sh
rev=`git rev-parse HEAD`
echo "changed_files=\`git diff --name-only $rev HEAD\`
touch --no-create \`echo \$changed_files\`" >> ${IEMIC_BUILD_DIR}/restore_timestamps.sh
