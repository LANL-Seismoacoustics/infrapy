#!/bin/bash
echo "Checking github for updates..."
git pull
echo 
echo "Checking for updates in the yml file"
conda env update --name infrapy_env --file infrapy_env.yml --prune
