#/bin/bash

set -e

## run on wenjie's droplets
build_dir=$(pwd)

cd $HOME/wenjie/wenjie-stat.me/
git pull origin master
cp -r $build_dir/docs/* $HOME/wenjie/wenjie-stat.me/static/reda/
git add -u docs/
git commit -m "deploy $CI_COMMIT_SHORT_SHA by gitlab-runner"
git push origin master
