#!/usr/bin/env bash
ffmpeg -r 2 -f image2 -i pca_gif%02d.png -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -vcodec libx264 -crf 5 -pix_fmt yuv420p test.mp4
