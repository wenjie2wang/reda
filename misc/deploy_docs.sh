#/bin/bash

set -e

## run on wenjie's droplets
pkg=$(grep "Package" DESCRIPTION | awk '{print $NF}')
build_dir=$(pwd)
target_dir=$HOME/wenjie/wenjie-stat.me/static/$pkg
tmp_log=.git_status.log
cd $HOME/wenjie/wenjie-stat.me/
git checkout -f
git checkout master
git pull origin master
mkdir -p $target_dir
cp -r $build_dir/docs/* $target_dir
git status > $tmp_log
if egrep -q "modified:[ ]+static/$pkg/" $tmp_log
then
    git add -u static/$pkg/
    git commit -m "deploy $pkg $CI_COMMIT_SHORT_SHA by gitlab-runner"
    git push origin master
else
    printf "The docs was not updated.\n"
fi
rm $tmp_log
