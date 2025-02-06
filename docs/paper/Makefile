SECTIONS = abstract ack body conclusion future_work intro \
	performance related titlepage
FIGURES = figures/arbfn_protocol.png figures/ffield_protocol.png

main.pdf:	main.tex refs.bib $(SECTIONS:%=sections/%.tex) $(FIGURES)
	pdflatex -shell-escape $<
	# bibtex main.aux
	pdflatex -shell-escape $<
	pdflatex -shell-escape $<

figures/%.png:	figures/%.plantuml
	plantuml $<

.PHONY:	clean
clean:
	find . \( -iname "*.aux" -or -iname "*.pdf"  \) \
		-exec rm -f "{}" \;
	rm -f *.aux *.bbl *.blg *.listing *.log *.out *.pdf
