all: figures

figures: faintpm_bd_qso.ps faintpm_test_plot.ps
	\cp -fv faintpm*.ps ../../tex/faintpm/
	touch ./figures

faintpm_bd_qso.ps: faintpm_bd.fits faintpm_highz.fits *.pro
	echo faintpm_bd_qso_plot | idl

faintpm_test_plot.ps: faintpm_test.fits *.pro
	echo faintpm_test_plot | idl

faintpm_bd.fits:
	echo "faintpm_highz_qsos,/bd" | idl

faintpm_highz.fits:
	echo "faintpm_highz_qsos" | idl
