#!/bin/bash
_now=$(date +'%H_%M_%S')
ffmpeg -f image2 -i %d.png video/video-$_now.avi -y
ffmpeg -i video/video-$_now.avi -pix_fmt rgb24 video/out_final-$_now.gif -y
find . -name '*.jpg' -delete
