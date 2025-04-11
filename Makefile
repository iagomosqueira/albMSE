
report/presentation_tcmp.pdf: presentation_tcmp.md
	R -e "rmarkdown::render('presentation_tcmp.md', output_dir='report')"

report/presentation_mse.pdf: presentation_mse.md
	R -e "rmarkdown::render('presentation_mse.md', output_dir='report')"

report/presentation_abc.pdf: presentation_abc.tex
	export TEXINPUTS=.:./boot/initial/report/:
	pdflatex -output-directory report presentation_abc.tex

report/IOTC-2025-TCMP09-04_Mosqueira-Hillary.pdf: report_tcmp.Rmd
	R -e "rmarkdown::render('report_tcmp.Rmd', output_dir='report', output_file='IOTC-2025-TCMP09-04_Mosqueira-Hillary.pdf')"

report/IOTC-2025-TCMP09-04_Mosqueira-Hillary-FR.pdf: report_tcmp_FR.Rmd
	R -e "rmarkdown::render('report_tcmp_FR.Rmd', output_dir='report', output_file='IOTC-2025-TCMP09-04_Mosqueira-Hillary-FR.pdf')"
