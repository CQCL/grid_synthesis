#!/bin/sh

eval "$(ssh-agent)"
ssh-add ~/.ssh/github_nihargargava_solovay
ssh-ad ~/.ssh/nihargargava_github_minkowski



git pull
git add .
git commit -m "Nihar has commited"
git push
