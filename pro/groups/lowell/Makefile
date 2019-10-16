
all: illustrate_nonparam_1.jpg illustrate_nonparam_2.jpg illustrate_nonparam_3.jpg illustrate_nonparam_4.jpg

illustrate_nonparam_%.ps: illustrate_nonparam.pro
	echo "illustrate_nonparam,plotfilename='illustrate_nonparam'" | idl

%.jpg: %.ps
	convert $< $@
