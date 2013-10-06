# Tell 'make' not to look for files with the following names
.PHONY : all

# By default, just build the docs
all :
	doxygen
	doxygen2docstring.py xml/index.xml SphericalFunctions_Doc.i
	perl -p -i -e 's/(%feature\("docstring"\) .*)((?:<(?:[^<>]++|(?1))*>))/$$1 . join("",split(" ",$$2))/ge; s/%feature\("docstring"\).* (\S+) """/%feature\("docstring"\) $$1 """/' SphericalFunctions_Doc.i
	perl -p -i -e 's/SphericalFunctions::SphericalFunctions::/SphericalFunctions::/g' SphericalFunctions_Doc.i
