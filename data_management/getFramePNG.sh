#!/bin/sh

echo "What is the base name of the file to process (no extension)?"
read infile
echo "What is the extension (e.g. .mp4)?"
read ext
echo "What time should image be extracted?"
read time

inname="${infile}${ext}"
outname="${infile}.png"

# Get frame from video
ffmpeg -ss ${time} -i ${inname} -s 1000x1000 -vframes 1 -f image2 ${outname}
