#!/bin/bash
_now=$(date +'%H_%M_%S')
ffmpeg -f image2 -i %d.png video/video-$_now.avi -y
ffmpeg -i Video_Handled/video-$_now.avi -pix_fmt rgb24 Video_Handled/out_final-$_now.gif -y
find . -name '*.jpg' -delete