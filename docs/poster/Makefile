# `make`-able figures from `../papers/`
FIGS = figures/arbfn_protocol.png figures/ffield_protocol.png

outline.pdf:	outline.md
	$(MAKE) -C ../paper $(FIGS)
	pandoc $< -o $@
