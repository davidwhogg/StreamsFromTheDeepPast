all: Methods.pdf Derivatives.pdf Fitting.pdf

Methods.pdf: Methods.tex
	pdflatex Methods.tex
	pdflatex Methods.tex
	pdflatex Methods.tex

Fitting.pdf: Fitting.tex
	pdflatex Fitting
	pdflatex Fitting
	pdflatex Fitting

Derivatives.pdf: Derivatives.tex
	pdflatex Derivatives
	pdflatex Derivatives
	pdflatex Derivatives

EllipseTest.png:
	python ../../simulation/Simulator.py
clean:
	rm *.aux *.log *.gz\(busy\)