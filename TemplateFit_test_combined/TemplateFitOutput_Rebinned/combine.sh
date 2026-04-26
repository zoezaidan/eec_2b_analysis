i=0
paste <(ls prefits/*.png | sort -V) <(ls postfits/*.png | sort -V) | \
while read f1 f2; do
  printf -v num "%03d" $i

  montage "$f1" "$f2" \
    -tile 2x1 \
    -geometry 1200x1200+0+0 \
    "combined/page_$num.png"

  ((i++))
done
