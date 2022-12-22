#!/bin/sh

eval "$(ssh-agent)"
ssh-add ~/.ssh/github_nihargargava_solovay
ssh-add ~/.ssh/github_gargava
ssh-add ~/.ssh/nihargargava_github_minkowski



git pull
git add .
git commit -m "Nihar has commited"
git push
