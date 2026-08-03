# Build the methods document.  Figures are pulled from ../out/ , so run
#   python -m fit.plot_fit --sample B
# first if they are missing or stale.
YbYAG_fit_method.pdf: method.tex ../out/fig_B_nk.png ../out/fig_B_yb_band.png ../out/fig_B_fit_vs_data.png
	pdflatex -interaction=nonstopmode method.tex >/dev/null
	pdflatex -interaction=nonstopmode method.tex >/dev/null
	mv method.pdf $@
	rm -f method.aux method.log method.out method.toc

clean:
	rm -f method.aux method.log method.out method.toc method.pdf

.PHONY: clean
