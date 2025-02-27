GIFS = 3_no_wall_effect.gif 4_wall_effect_5.gif \
	5_wall_effect_10.gif 6_wall_effect_20.gif \
	7_efficient_wall_effect_15.gif \
	8_efficient_wall_effect_90.gif topology.png

pres.pdf:	pres.md $(GIFS)
	pandoc -t beamer $< -o $@

%.gif:	%.avi
	ffmpeg -i $< -t 10 $@

%.png:	%.gv
	dot -Tpng -o $@ $<

clean:
	rm -f $(GIFS) pres.pptx
