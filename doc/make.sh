pandoc --template=template.html --mathjax --toc --smart --metadata=title:milonga ../README.md -o html/index.html
pandoc --template=template.html --mathjax --toc --smart ../INSTALL.md -o html/install.html
# pandoc --template=template.html --mathjax --toc --smart FAQ.md -o html/FAQ.html
# TODO citations
# pandoc --template=template.html --mathjax --toc --smart --number-sections description.md -o html/description.html
# pandoc --template=template.html --mathjax --toc --smart --number-sections realbook.md -o html/realbook.html

# pandoc --template=template.texi wasora.md -o wasora.texi


# cd ../examples
# ./showcase-html.sh > ../doc/showcase.md
# cd ../doc
# 
# pandoc --template=template.html --mathjax --toc --smart showcase.md -o html/showcase.html
# cp ../examples/lorenz.svg html

