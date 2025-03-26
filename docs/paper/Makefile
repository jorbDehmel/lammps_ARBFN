FIGURES = figures/arbfn_protocol.png \
	figures/ffield_protocol.png figures/ffield_every_protocol.png

main.pdf:	main.tex refs.bib $(FIGURES)
	pdflatex $<
	bibtex main.aux
	pdflatex $<
	pdflatex $<

figures/%.png:	figures/%.plantuml
	plantuml $<

.PHONY:	clean
clean:
	find . \( -iname "*.aux" -or -iname "*.pdf"  \) \
		-exec rm -f "{}" \;
	rm -f *.aux *.bbl *.blg *.listing *.log *.out *.pdf
