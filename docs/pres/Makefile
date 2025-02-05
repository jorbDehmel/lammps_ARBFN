GIFS = 0_no_fix.gif 1_dampened.gif 2_halted.gif \
	3_no_wall_effect.gif 4_wall_effect_5.gif \
	5_wall_effect_10.gif 6_wall_effect_20.gif \
	7_efficient_wall_effect_15.gif \
	8_efficient_wall_effect_90.gif

pres.pdf:	pres.md $(GIFS)
	pandoc -t beamer $< -o $@

%.gif:	%.avi
	ffmpeg -i $< -t 10 $@

clean:
	rm -f $(GIFS) pres.pptx
