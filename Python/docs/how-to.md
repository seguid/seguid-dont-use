### How to generate these docs


Make an empty doc folder, then:

	sphinx-quickstart

This generates lots of files including

	conf.py

Add this to conf.py:

	extensions = ['sphinx.ext.autodoc',
		          'sphinx.ext.coverage',
		          'sphinx.ext.napoleon']


Then:

	sphinx-apidoc --force --no-toc --no-headings --output-dir . ..




sphinx-build . _build/html&&xdg-open _build/html/index.html
