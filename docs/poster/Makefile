FIGS = figures/arbfn_protocol.png figures/ffield_protocol.png

outline.pdf:	outline.md
	$(MAKE) -C ../arxiv $(FIGS)
	pandoc $< -o $@
