#!/usr/bin/env bash
ffmpeg -r 2 -f image2 -s 500x500 -i spatial_gif%02d.png -vcodec libx264 -crf 5 -pix_fmt yuv420p test.mp4
