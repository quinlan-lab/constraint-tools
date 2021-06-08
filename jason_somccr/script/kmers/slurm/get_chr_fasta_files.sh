#!/bin/bash

url="$1"
filename="$2"

echo "$url"
echo "$filename"

curl "$url" | gunzip -c > "$filename"

samtools faidx "$filename"