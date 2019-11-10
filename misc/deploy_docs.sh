#/bin/bash

set -e

## run on wenjie's droplets
build_dir=$(pwd)
cd $HOME/wenjie/wenjie-stat.me/
git checkout -f
git checkout master
git pull origin master
cp -r $build_dir/docs/* $HOME/wenjie/wenjie-stat.me/static/reda/
git add -u static/reda/
git commit -m "deploy $CI_COMMIT_SHORT_SHA by gitlab-runner"
git push origin master
