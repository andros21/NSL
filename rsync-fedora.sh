#! /usr/bin/env bash
rsync -aRvh --progress --stats --exclude='*.x' --exclude='*.o' --exclude='*.gch' ./TSP ./RC ./10 fedora:/home/fedora/NSL 
