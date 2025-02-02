GIFS = 0_no_fix.gif 1_dampened.gif 2_halted.gif \
	3_no_wall_effect.gif 4_wall_effect_5.gif \
	5_wall_effect_10.gif 6_wall_effect_20.gif

pres.pdf:	pres.md $(GIFS)
	pandoc -t beamer $< -o $@

%.gif:	proof_of_concept/%.avi
	ffmpeg -i $< -t 10 $@

clean:
	rm -f $(GIFS) pres.pptx
