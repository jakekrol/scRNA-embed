gargs_linux --log g.log -o "Rscript --no-init-file rds2mtx.R {0} {1} && sleep 3 && ./mtx2h5ad.py {2} {3} {4} {5}" < rds2mtx2h5ad.input
