TARGET	=  main
NAME	=  nonlinear-control_miniproject
MODE	=  -interaction=batchmode --jobname=${NAME}
GARBAGE =  *.aux *.bbl *.blg *.log *.lox *.toc *.out *.lof *.lot *.nav *.snm main.pdf *.acn *.acr *.alg *.glg *.glo *.gls *.ist *.backup *.dvi rapport.pdf *.tdo network-note.pdf

# Colour definitions
CGREEN	=  "\033[1;92m"
CNORMAL	=  "\033[0;37m"
CYELLOW =  "\033[1;33m"

.PHONY: pdf clean show

all: pdf clean

pdf:
	@echo ${CGREEN} "Building miniproject" ${CNORMAL}
	@echo ${CYELLOW} "make pdf" ${CNORMAL}
	pdflatex ${MODE} ${TARGET}.tex
	pdflatex ${MODE} ${TARGET}.tex
	@echo ${CGREEN} "Build done" ${CNORMAL}
	@echo ""

clean:
	@echo ${CYELLOW} "Cleaning up garbage files" ${CNORMAL}
	@find . | grep ".backup" | xargs rm -rf
	@rm -f ${GARBAGE}
	@echo ${CGREEN} "Done" ${CNORMAL}

show:
	xdg-open ${NAME}.pdf
