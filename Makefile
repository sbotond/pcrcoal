
.PHONY: ct cat com push gt checkpkg clean remove aareload
PKG=pcrcoal_1.1.tar.gz

ct:
	git log --graph
cat: *.R
	(rm PCRcoalSource.R;cat *.R > PCRcoalSource.R;true)
com: *.R
	git commit -a
push:
	git push --all
fetch:
	git fetch --all
gt:
	gitk --all
pkg: *.R cat
	(rm pkg/R/*.R;true)
	(rm PCRcoalSource.R;true)
	cp *.R pkg/R/
	R CMD build pkg
checkpkg: pkg 
	R CMD check --as-cran $(PKG)
clean:
	(rm *.log; rm $(PKG);rm -r ./pcrcoal.Rcheck; rm ./pkg/R/*.R;true ) 2>&1 > /dev/null
inst: pkg
	R CMD INSTALL	$(PKG)
remove:
	R CMD REMOVE pcrcoal

