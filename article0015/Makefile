.PHONY: default all
default: all
all: article0015

PYCULGI = python-culgi-6.0.0-kml

.PHONY: article0015
article0015: article0015-fig2.png article0015-fig3.png article0015-fig4-single.png article0015-fig4.png article0015-fig5.png article0015-fig6.png article0015-fig7.png
article0015-fig4-single.png: article0015.py
	"$(PYCULGI)" article0015.py fig4 single
article0015-%.png: article0015.py
	"$(PYCULGI)" article0015.py $*